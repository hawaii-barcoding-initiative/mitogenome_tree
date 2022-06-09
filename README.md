# mitogenome_tree
Mitogenome phylogenetic tree workflow. This does two things

1. `tree_prep.R` creates the concatenated alignment and partitionfinder config files (but it won't run partitionfinder or build the trees for you) 
2. Once you've run the MrBayes & RAXML analyses, come back here with your consensus trees (in data/tree/bayes.tre and data/tree/ml.tre) and run `plot_tree.R` to generate the tree figure used in the paper (give or take some minor tweaks in Illustrator)

### how to run this thing:

* install dependences
  open R and run `renv::restore()`
  intstall mafft however you want, e.g. `conda install mafft`
  make sure you've got python 3 with BioPython installed, like `pip3 install BioPython`

  you'll also need an anaconda instance with python 2 to run PartitionFinder, which is easily done:
  ```bash
  $ curl -LO https://repo.anaconda.com/archive/Anaconda2-5.3.1-Linux-x86_64.sh
  $ bash Anaconda2-5.3.1-Linux-x86_64.sh -b -p ~/anaconda
  $ source ~/anaconda/bin/activate
  # do stuff
  $ source ~/anaconda/bin/deactivate
  ```

* Replace the contents of data/api_key.txt with a valid NCBI entrez API key
