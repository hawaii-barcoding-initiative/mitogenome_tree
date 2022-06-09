library(tidyverse)
library(ggtree)
library(phytools)
library(here)
library(tidytree)


# load the data we already have about the species in our alignment
taxonomy <- read_csv(here("data","alignment_data.csv")) %>%
  mutate(label = str_replace(str_glue("{accession}-{species}")," ","_"))

# this is a manually-generated list of tree tip pairs that corresponds to family-level splits 
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
      bayes = treeio::read.mrbayes(here("data","tree",str_glue("{.x}.tre"))),
      ml = treeio::read.raxml(here("data","tree",str_glue("{.x}.tre")))
    )
    # reroot the tree on the flat sharks
    # first, find the parent node so we know where to reroot
    rn <- MRCA(tree,c("NC_024102-Gymnura_poecilura", "OK104094-Gymnura_altavela"))
    # calculate the new edge length
    len <- tree@phylo$edge.length[tree@phylo$edge[, 2] == rn]
    # reroot the tree with the new edge length
    tree@phylo <- reorder(reroot(tree@phylo,rn,len/2))
    # get 
    supp <- fam_combos %>%
      map_dfr(~{
        node <- MRCA(tree,c(.x[1],.x[2]))
        list(node1=.x[1],node2=.x[2],parent=node,level="family")
      }) 
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
                  list(node1=.x[1],node2=.x[2],parent=MRCA(trees$bayes$tree,c(.x[1],.x[2])),level="genus")
                }) %>%
                distinct(parent,.keep_all = TRUE)
            }
          }) %>%
          map_dfr(~.x)
      ) 
    fam <- taxonomy %>%
      group_by(family) %>% 
      summarise(node = MRCA(tree,label))
    list(tree=tree,support_nodes=supp,families=fam)
  })


  
node_support <- trees$bayes$support_nodes %>%
  inner_join(trees$ml$support_nodes,by = c("node1" = "node1","node2" = "node2"),suffix=c("_bayes","_ml")) %>%
  inner_join(as_tibble(trees$bayes$tree) %>% select(node,prob),by=c("parent_bayes" = "node")) %>%
  inner_join(as_tibble(trees$ml$tree) %>% select(node,bootstrap),by=c("parent_ml" = "node")) %>%
  mutate(across(c(prob,bootstrap),as.numeric)) %>%
  select(node=parent_bayes,bootstrap=bootstrap,level=level_bayes)

td <- taxonomy %>%
  separate(species,into=c("genus","species"),sep=" ") %>%
  mutate(abbrev = str_c(str_sub(genus,1,1),". ",species)) %>%
  unite("species",genus,species,sep=" ",remove=FALSE) %>%
  mutate(color = case_when(
    in_study ~ "Sequenced in this study",
    TRUE ~ "Downloaded from MiFish"
  )) %>%
  select(label,species,abbrev,color)

quartz()
bayes_tree <- trees$bayes$tree %>%
  as_tibble() %>%
  left_join(node_support,by="node") %>%
  mutate(prob=as.numeric(prob)) %>%
  mutate(support = str_c(prob_percent,"/",bootstrap)) %>%
  # mutate(support = str_pad(support, width=max(str_length(support),na.rm=TRUE), side="left")) %>%
  as.treedata()
(treep <- ggtree(bayes_tree) %<+% td +
  geom_nodelab(
    aes(label=support,subset = node %in% node_support$node & level == "family"),
    nudge_x = -0.11,
    nudge_y = 0.5
  ) +
  geom_nodelab(
    # geom="label",
    nudge_x = -0.11,
    nudge_y = 0.5,
    aes(label=support,subset = node %in% node_support$node & level == "genus"),
    size=2
  ) +
  geom_tiplab(aes(label=species,color=color),fontface="italic") + 
  # geom_cladelab(mapping=aes(node=node,label=family),data=trees$bayes$families,align=TRUE,hjust=0,offset=0.27) + 
  geom_cladelab(mapping=aes(node=node,label=family),barsize=2,data=trees$bayes$families,hjust=0,offset=0.4,fontface="bold",align=TRUE) +
  scale_color_manual(values=c("black","firebrick")) +
  xlim(-0.02,4) +
  theme(legend.position = "none"))
ggsave("data/tree/tree.pdf",device="pdf",height=8.5,width=11,units="in")
