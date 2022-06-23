# mitogenome_tree


Mitogenome phylogenetic tree workflow. This does two things

1. `tree_prep.R` creates the concatenated alignment and partitionfinder config files (but it won't run partitionfinder or build the trees for you) 
2. Once you've run the MrBayes & RAXML analyses, come back here with your consensus trees (in data/output/bayes.tre and data/output/ml.tre) and run `plot_tree.R` to generate the tree figure used in the paper.

### how to run this thing:

* install dependencies  
  open R and run `renv::restore()`  
  install mafft however you want, e.g. `conda install mafft`  
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

* Run the code in `tree_prep.R`

* Use the output file `partition_finder.cfg` to run PartitionFinder 2

  ```bash
  user@computer ~/code/mitogenome_tree $ cd data/partition
  user@computer ~/code/mitogenome_tree $ source ~/anaconda/bin/activate
  user@computer ~/code/mitogenome_tree $ python ~/path/to/PartitionFinder.py .
  .... a bunch of stuff happens ....
  ```

* Use data/partitions/analysis/best_models.txt and the output file `concatenated_alignment.nex` to make your MrBayes nexus file 

* Use data/partitions/analysis/best_models.txt and the output file `concatenated_alignment.phy` to run RAxML. 

* Run MrBayes and RAxML

* Copy the output trees to data/output/bayes.tre and data/output/ml.tre and run `plot_tree.R`

* Now you have the tree plot from the paper
