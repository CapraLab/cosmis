# COSMIS
COSMIS is a novel framework for quantifying the 3D spatial constraint of amino acid sites in the human proteome. If you find COSMIS useful in your work, please consider citing the COSMIS paper: 
* Li, B., Roden, D.M., and Capra, J.A. (2021). The 3D spatial constraint on 6.1 million amino acid sites in the human proteome. bioRxiv. doi: https://doi.org/10.1101/2021.09.15.460390

## What's in each folder?

### cosmis
The `cosmis` folder contains utility code that the top level cosmis code `cosmis.py`, `cosmis_batch.py`, and `cosmis_sp.py` depend on.

### fig_scripts
The `fig_scripts` folder contains standlone R code that can be run to reproduce figures in the main text and supplementary document of the manuscript.

### data
The `data` folder constains all raw and processed data that code in the `fig_scripts` will load to reproduce the figures.

### scores
The `scores` folder constains precomputed scores for all 16,533 proteins of the human reference proteome currently covered by the framework.

### scripts
The `scripts` folder contains all scripts written in Python that were called to obtain processed datasets in the `data` folder.

## More descriptions and tutorials are under construction
