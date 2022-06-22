# TODO: ML support values are correct but bayesian ones might be wonky
library(tidyverse)
library(ggtree)
library(phytools)
library(here)
library(tidytree)
library(svglite)
library(ochRe)
library(treeio)


source("reroot.R")


# load the data we already have about the species in our alignment
taxonomy <- read_csv(here("data","alignment_data.csv")) %>%
  mutate(label = str_replace(str_glue("{accession}-{species}")," ","_"))

# this is a manually-generated list of tree tip pairs that corresponds to
# family-level splits in both trees
fam_combos <- list(
  c("NC_004417-Gymnothorax_kidako","NC_023980-Saurida_microlepis"),
  c("NC_004417-Gymnothorax_kidako","NC_024102-Gymnura_poecilura"),
  c("NC_023980-Saurida_microlepis","NC_011180-Ablennes_hians"),
  c("NC_011180-Ablennes_hians","NC_026718-Decapterus_macarellus"),
  c("NC_011180-Ablennes_hians","NC_004379-Lota_lota"),
  c("NC_011180-Ablennes_hians", "NC_003189-Myripristis_berndti"),
  c("NC_026718-Decapterus_macarellus","NC_009870-Chaetodon_auripes"),
  c("NC_009870-Chaetodon_auripes","NC_009865-Ostracion_immaculatus"),
  c("NC_009865-Ostracion_immaculatus","NC_010978-Canthigaster_coronata")
)

  # load the two trees and reroot them on the flat sharks, return them in a list
trees <- c("bayes","ml") %>%
  set_names %>%
  purrr::map(~{
    tree <- switch(
      .x, # figure out how to load the tree files
      bayes = {
        treeio::read.mrbayes(here("data","output",str_glue("{.x}.tre"))) %>%
          as_tibble() %>%
          mutate(node = as.numeric(node)) %>%
          as.treedata()
        # t@data %>%
        #   mutate(node = as.numeric(node))
      },
      ml = treeio::read.raxml(here("data","output",str_glue("{.x}.tre")))
    )
    # reroot the tree on the flat sharks
    # first, find the parent node so we know where to reroot
    rn <- MRCA(tree,c("NC_024102-Gymnura_poecilura", "OK104094-Gymnura_altavela"))
    # calculate the new edge length
    len <- tree@phylo$edge.length[tree@phylo$edge[, 2] == rn]
    # reroot the tree with the new edge length
    # we use a custom function that keeps all the associated node support values correct
    tree <- reroot_treedata(tree,node=rn,len=len)
    # get the nodes at family-level splits and make it into a data frame. this is
    # a bunch of fairly convoluted stuff I did when I originally planned to show
    # support at only family- and genus-level nodes, but it contains logic that
    # allows me to combine the support for the two trees, so I'm leaving it in
    # rather than trying to figure out how to remove it without breaking
    # everything downstream
    supp <- fam_combos %>%
      map_dfr(~{
        node <- MRCA(tree,c(.x[1],.x[2]))
        list(node1=.x[1],node2=.x[2],parent=node,level="family")
      }) 
      # append the genus-level splits
    supp <- supp %>%
      bind_rows(
        taxonomy %>%
          group_by(family) %>%
          group_map(~{
            gen_combos <- .x %>%
              distinct(genus,.keep_all = TRUE) %>%
              pull(label)
            if (length(gen_combos) > 1) {
              gen_combos %>%
                combn(2,simplify = FALSE) %>%
                map_dfr(~{
                  list(node1=.x[1],node2=.x[2],parent=MRCA(tree,c(.x[1],.x[2])),level="genus")
                }) %>%
                distinct(parent,.keep_all = TRUE)
            }
          }) %>%
          map_dfr(~.x)
      )
    
    # this allows us to label clades by family
    fam <- taxonomy %>%
      group_by(family) %>% 
      summarise(node = MRCA(tree,label))
    list(tree=tree,support_nodes=supp,families=fam)
  })

  
# generate a combined data frame of node support values for both trees
node_support <- trees$bayes$support_nodes %>%
  inner_join(trees$ml$support_nodes,by = c("node1" = "node1","node2" = "node2"),suffix=c("_bayes","_ml")) %>%
  inner_join(as_tibble(trees$bayes$tree) %>% select(node,prob),by=c("parent_bayes" = "node")) %>%
  inner_join(as_tibble(trees$ml$tree) %>% select(node,bootstrap),by=c("parent_ml" = "node")) %>%
  mutate(across(c(prob,bootstrap),as.numeric)) %>%
  select(node=parent_bayes,bootstrap=bootstrap,level=level_bayes)

# this is a bit of extra data to associate with the tree when plotting it
# some of it isn't used but I left the code in
# mostly what it does it make the species in the study bold & italic
td <- taxonomy %>%
  separate(species,into=c("genus","species"),sep=" ") %>%
  mutate(abbrev = str_c(str_sub(genus,1,1),". ",species)) %>%
  unite("species",genus,species,sep=" ",remove=FALSE) %>%
  mutate(face = case_when(
    in_study ~ "bold.italic",
    TRUE ~ "italic"
  )) %>%
  select(label,species,abbrev,in_study,face) 

# names of factor levels to split node support into
prob_labels <- c("<60%","60–75%","75–95%",">95%")

# let's rebuild a unified tree with all the data we need
bayes_tree <- trees$bayes$tree %>%
  as_tibble() %>%
  left_join(node_support,by="node") %>% # pull in bootstrap support from the other tree
  mutate(prob=as.numeric(prob)) %>%     # make prob numeric
  mutate(
    posterior_bins = cut(               # cut posterior into factor levels
      prob,
      c(-Inf,0.60,0.75,0.95,Inf), 
      labels=prob_labels
    ),
    boot_bins = cut(                    # cut bootstrap into factor levels
      bootstrap,
      c(-Inf,60,75,95,Inf),
      labels=prob_labels
    )
  ) %>%
  as.treedata()

fontsize <- 3

# plot the thing
(treep <- ggtree(bayes_tree) %<+% td +
  # make node symbols for bootstrap values <= 95 at family- and genus-level splits
  geom_nodepoint(aes(fill=boot_bins,subset=bootstrap <= 95),shape=23,alpha=1,size=3) +
  # label posterior probabilities <= 95% for family- and genus-level splits
  # there are two calls because we nudge the label over less if there's a symbol taking up space there
  geom_nodelab(aes(label=prob_percent,subset = prob <= 0.95 & (bootstrap > 95 | is.na(bootstrap)) & (level %in% c("family","genus"))  ),hjust = -0.5,nudge_x = -0.03,size=fontsize) +
  geom_nodelab(aes(label=prob_percent,subset = prob <= 0.95 & bootstrap <= 95 & (level %in% c("family","genus"))),hjust = -0.5,nudge_x = -0.01,size=fontsize) +
  geom_treescale(x=2) +  # show a branch length scale
  geom_tiplab(aes(label=species,color=in_study,fontface=face)) +  # label species
  # label families
  geom_cladelab(mapping=aes(node=node,label=family),barsize=2,data=trees$bayes$families,hjust=0,offset=0.6,fontface="bold",align=TRUE) +
  scale_color_manual(values=c("black","firebrick")) + # color species labels
  scale_fill_ochre(palette="emu_woman_paired",name="ML bootstrap support",reverse=TRUE) +  # color node symbols
  guides(color="none") + # don't show a legend for the font color
  xlim(-0.02,3.6) +
  labs(fill="ML bootstrap support") +
  theme(legend.position = c(0.1,0.85)))

# save in two different formats, cairo_pdf lets us use unicode characters (there's an en-dash in the legend)
ggsave(plot=treep,here("data","output","Figure 3.svg"),device="svg",height=8.5,width=11,units="in")
ggsave(plot=treep,here("data","output","Figure 3.pdf"),device=cairo_pdf,height=8.5,width=11,units="in")
