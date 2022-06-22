new_labels <- function (tree) 
{
  node2label_old <- tree %>% as_tibble() %>% dplyr::select(c("node", 
                                                             "label"))
  if (inherits(tree, "treedata")) {
    tree <- tree@phylo
  }
  tree$tip.label <- paste0("t", seq_len(Ntip(tree)))
  tree$node.label <- paste0("n", seq_len(Nnode(tree)))
  node2label_new <- tree %>% as_tibble() %>% dplyr::select(c("node", 
                                                             "label"))
  old_and_new <- node2label_old %>% dplyr::inner_join(node2label_new, 
                                                      by = "node") %>% dplyr::rename(old = "label.x", new = "label.y")
  return(list(tree = tree, node2old_new_lab = old_and_new))
}

map_nodes <- function (oldtree, newtree) 
{
  treelab1 <- oldtree %>% as_tibble() %>% dplyr::select(c("node", 
                                                          "label"))
  treelab2 <- newtree %>% as_tibble() %>% dplyr::select(c("node", 
                                                          "label"))
  node_map <- dplyr::inner_join(treelab1, treelab2, by = "label") %>% 
    dplyr::select(c("node.x", "node.y")) %>% dplyr::rename(c(old = "node.x", 
                                                             new = "node.y"))
  return(node_map)
}

new_tree <- function (tree, node2old_new_lab) 
{
  treeda <- tree %>% as_tibble()
  treeda1 <- treeda %>% dplyr::filter(.data$label %in% node2old_new_lab$new)
  treeda2 <- treeda %>% dplyr::filter(!(.data$label %in% node2old_new_lab$new))
  treeda1$label <- node2old_new_lab[match(treeda1$label, node2old_new_lab$new), 
                                    "old"] %>% unlist(use.names = FALSE)
  treeda <- rbind(treeda1, treeda2)
  tree <- treeda[order(treeda$node), ] %>% as.phylo()
  return(tree)
}

reroot_treedata <- function (phy, node = NULL, len)
{
  # if (!missing(outgroup) && is.character(outgroup)) {
  #   outgroup <- match(outgroup, phy@phylo$tip.label)
  # }
  # if (!edgelabel) {
  #   message("The use of this method may cause some node data to become incorrect (e.g. bootstrap values) if 'edgelabel' is FALSE.")
  # }
  res <- new_labels(tree = phy)
  tree <- res$tree
  node2oldnewlab <- res$node2old_new_lab
  re_tree <- reorder(reroot(tree,node,len/2))
  # re_tree <- ape::root(tree, outgroup = outgroup, node = node, 
  #                 edgelabel = edgelabel)
  node_map <- map_nodes(tree, re_tree)
  n.tips <- Ntip(re_tree)
  phy@phylo <- new_tree(tree = re_tree, node2old_new_lab = node2oldnewlab)
  update_data <- function(data, node_map) {
    cn <- colnames(data)
    cn <- cn[cn != "node"]
    data <- dplyr::inner_join(data, node_map, by = c(node = "old")) %>% 
      dplyr::select(c("new", cn)) %>% dplyr::rename(node = .data$new)
    root <- data$node == (n.tips + 1)
    data[root, ] <- NA
    data[root, "node"] <- n.tips + 1
    return(data)
  }
  if (nrow(phy@data) > 0) {
    phy@data <- update_data(phy@data, node_map)
  }
  if (nrow(phy@extraInfo) > 0) {
    phy@extraInfo <- update_data(phy@extraInfo, node_map)
  }
  return(phy)
}

