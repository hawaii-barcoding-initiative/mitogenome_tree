library(tidyverse)
library(here)
library(fs)
library(rentrez)
library(parallel)
library(httr)
library(XML)
library(ape)

# TODO: handle rRNAs

# NCBI entrez API key
api_key <- "<api_key>"

# load study species
study_species <- read_csv(here("data","study-species.csv")) %>%
  mutate(species = str_glue("{genus} {species}")) %>%
  select(order,family,genus,species,accession)
# load mitofish species

mitofish_url <- "http://mitofish.aori.u-tokyo.ac.jp/species/all.html"
mt_species <- readHTMLTable(mitofish_url,header=T,as.data.frame = T, which=1,trim = T) %>% 
  as_tibble() %>%
  unite("Species",Genus,Species,sep = " ",remove=FALSE) %>%
  mutate(accession = NA_character_) %>%
  # select(-sci_name,sci_name) %>%
  rename_with(~tolower(.x)) %>% 
  # they're missing a family specifically for Gymnura for some reason, so we gotta fix that
  mutate(
    family = case_when(
      genus == "Gymnura" ~ "Gymnuridae",
      TRUE ~ family
    )
  ) %>%
  select(order,family,genus,species,accession)

# mt_species <- read_csv(here("data","mitofish-species.csv")) %>%
#   mutate(sci_name = str_glue("{genus} {species}"))

# get taxon information from rfishbase
# all_taxa <- load_taxa() %>%
#   filter(Species %in% unique(c(study_species$sci_name,mt_species$sci_name))) %>%
#   collect()

# download the mifish sequences. they don't provide accession numbers in their data table
# or else I could avoid doing this. it's a lot to download just for a few numbers
download <- TRUE
if (download) {
  mitofish <- "http://mitofish.aori.u-tokyo.ac.jp/files/mitogenomes.zip"
  mf_dir <- dir_create(here("data","mitogenomes"))
  mf_file <- path(mf_dir,"mitogenomes.zip")
  cat("downloading mitofish sequences (unfortunately, there doesn't seem to be a better way to get their accession numbers)\n")
  GET(mitofish,write_disk(mf_file,overwrite=TRUE),progress())
  utils::unzip(mf_file,exdir = mf_dir,overwrite = TRUE)  
}


# get lineage data for mitofish species
# mt_species <- mt_species %>%
#   inner_join(all_taxa,by=c("sci_name" = "Species")) %>%
#   select(superclass=SuperClass,class=Class,order=order,family=Family,genus=Genus,species=sci_name) %>%
#   mutate(mitogenome=TRUE)

# get lineage data for study species and count how many of each level are in mitofish
sp <- study_species %>%
  # inner_join(all_taxa,by=c("sci_name" = "Species")) %>%
  # select(superclass=SuperClass,class=Class,order=order,family=Family,genus=Genus,species=sci_name) %>%
  mutate(
    across(
      order:species,
      function(taxon) {
        map_dbl(taxon,~sum(mt_species %>% pull(.y) %in% .x),cur_column())
      },
      .names="n_{.col}"
    )
  )

# from looking at the summary, we know that family is the level that represents everything in our study
# get list of species in this study, genera represented in mitofish, and families represented in mitofish
spp <- sp$species
genera <- sp %>% filter(n_genus > 0) %>% pull(genus)
families <- sp %>% filter(n_family > 0 & n_genus == 0) %>% pull(family)

# there are a ton of carangids so we're just gonna keep a few (presumably closely related) genera
# here we filter the list and try to find the filenames for all the fasta files
# we also extract the accession numbers from the filenames (we're gonna need those)
# you have to have downloaded all the mitofish mitogenomes. In fact, I think I'll do that up above
keep_carangids <- c("Seriola","Elagatis","Decapterus")
mitos <- mt_species %>%
  # select(order:species) %>%
  mutate(in_study=FALSE) %>%
  filter(species %in% spp | genus %in% genera | family %in% families) %>%
  filter(family != "Carangidae" | genus %in% keep_carangids) %>%
  bind_rows( sp %>% mutate(in_study=TRUE) ) %>%
  mutate(
    file_pattern = str_replace(species," ","_"),
    fasta = map_chr(file_pattern,~dir_ls(here("data","mitogenomes"),regexp=.x,recurse = TRUE)[1]),
    accession = case_when(
      # this is where we have to extract the accession numbers from the filename
      is.na(accession) ~ str_extract(path_file(fasta),"^((?:NC_)?[a-zA-Z0-9]+)"),
      TRUE ~ accession
    )
  ) %>%
  select(-starts_with("n_"),-fasta,-file_pattern)

# save the information we've collected so far
write_csv(mitos,here("data","alignment_data.csv"))

# let's pull down the full genbank entries fro all the accession numbers we've got
gb_file <- entrez_fetch(db="nuccore",id=mitos$accession,rettype="gb",api_key=api_key)
# save those into a flat file
write_lines(gb_file,here("data","tree_mitogenomes.gb"))

# now use this clever python script to separate out the protein-coding regions
system(str_glue("python3 {here('coding_sequences.py')} {here('data','tree_mitogenomes.gb')} > {here('data','mitogenome_cds.tab')}"))

# here's our alignment directories
pa_dir <- here("data","prealignment")
al_dir <- here("data","alignment")
dir_create(pa_dir)
dir_create(al_dir)

# this part could have just been done inside the python script but for
# whatever reason I'm doing it back in this R script
coding_regions <- read_tsv("data/mitogenome_cds.tab")

# now we create separate fastas for each coding region so we can align them independently
coding_regions %>%
  mutate(gene = str_replace(gene,"sythase","synthase"),len=str_length(nucleotides)) %>%
  group_by(gene) %>%
  group_walk(~{
    filename <- here("data","prealignment",str_glue("{.y$gene}.fasta"))
    contents <- .x %>%
      pmap_chr(~{
        row <- list(...)
        str_c(
          str_replace(str_glue(">{row$accession}-{row$species}")," ","_"),
          "\n",
          row$nucleotides,
          collapse = "\n"
        )
      })
    write_lines(contents,filename,append=FALSE)
  })

# now we can align the different coding regions (assumes you've got mafft installed)
cores <- detectCores()
# setup our command string for running mafft
cmd <- str_glue("for f in {path(pa_dir,'*.fasta')}; do mafft --thread {cores} --auto $f > {al_dir}/$(basename $f .fasta)_aligned.fasta; done")
# run mafft
system(cmd)

gene_types <- coding_regions %>%
  distinct(gene,type)

# now we want to concatenate them all and keep a record of their location
start <- 1
concatenated <- NULL
alignments <- dir_ls(al_dir,glob="*.fasta") %>%
  map2_dfr(seq_along(.),~{
    alignment <- read.dna(.x,"fasta")
    len <- dim(alignment)[2]
    gene <- str_extract(path_file(.x),"^(\\w+)(?=_aligned\\.fasta)")
    if (is.null(concatenated)) {
      concatenated <<- alignment
    } else {
      concatenated <<- cbind(concatenated,alignment,fill.with.gaps=TRUE)
    }
    row <- list("gene" = gene, "start" = as.numeric(start), "length" = as.numeric(len), "index" = .y)
    start <<- start + len
    return(row)
  })

# create partitionfinder.cfg file
cfg <- str_c(
  "### ALIGNMENT ###",
  str_glue("alignment = concatenated_alignment.phy;"),
  "",
  "## BRANCHLENGTHS: linked | unlinked ##",
  "branchlengths = linked;",
  "",
  "## MODELS OF EVOLUTION: all | allx | mrbayes | beast | gamma | gammai | <list> ##",
  "models = mrbayes;",
  "",
  "# MODEL SELECCTION: AIC | AICc | BIC #",
  "model_selection = aicc;",
  "",
  "## DATA BLOCKS: see manual for how to define ##",
  "[data_blocks]",
  str_c(pmap_chr(alignments,~{
    row <- list(...)
    gt <- gene_types %>%
      filter(gene == row$gene) %>%
      pull(type)
    positions <- ""
    if (gt[1] == "CDS") {
      positions <- map_chr(seq(3),function(i) {
        name <- str_glue("{row$gene}_codon{i}")
        start <- row$start + i- 1
        end <- row$start + row$length - 1
        str_glue("{name} = {start}-{end}\\3;")
      })
    } else if (gt == "rRNA") {
      name <- row$gene
      start <- row$start
      end <- row$start + row$length - 1
      positions <- str_glue("{name} = {start}-{end};")
    }
    str_c(positions,collapse="\n")
  }),collapse="\n"),
  "",
  "## SCHEMES, search: all | user | greedy | rcluster | rclusterf | kmeans ##",
  "[schemes]",
  "search = greedy;",
  sep="\n"
)

# write partitionfinder data
pf_dir <- here("data","partition")
dir_create(pf_dir)

# write concatenated alignment in phylip format
cc_file <- path(pf_dir,"concatenated_alignment.phy")
write.dna(concatenated,cc_file,format="sequential",nbcol=-1,colsep="")
write_lines(cfg,path(pf_dir,"partition_finder.cfg"))

