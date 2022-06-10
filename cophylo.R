library(dplyr)
library(ggtree)

# plot the two phylograms face to face
# you need to have run plot_tree.R first

x <- trees$bayes$tree
y <- trees$ml$tree
rot <- c(69,67,81,82,75,76,77,78,97,96,91,66,103)
p1 <- ggtree(x, ladderize=TRUE) 
p2 <- ggtree(y, ladderize=TRUE)
rot %>%
  purrr::walk(~{
    p2 <<- p2 %>% ggtree::rotate(.x)
  })

d1 <- p1$data
d2 <- p2$data

## reverse x-axis and 
## set offset to make the tree on the right-hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 1

pp <- p1 + geom_tree(data=d2) 
  # ggnewscale::new_scale_fill() 

dd <- bind_rows(d1, d2) %>% 
  filter(!is.na(label))

pp + geom_line(aes(x, y, group=label), data=dd, color='grey')
