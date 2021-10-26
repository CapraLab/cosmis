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

### structures
The `structures` folder contains all protein structure in the PDB format based on which COSMIS scores were computed.

### supplementary_tables
The `supplementray_tables` folder contains all supplementary tables referred to in the COSMIS paper.

## Using the COSMIS framework
It is recommended that interested users of the COSMIS framework download precomputed COSMIS scores from this repository. However, should you need to run COSMIS using custom-built protein structural models, or to compute COSMIS scores based on protein-protein complexes, please follow the following steps.

### Clone COSMIS
Clone COSMIS to a local directory.
```bash
git clone https://github.com/CapraLab/cosmis.git
```

### Set up conda environment
COSMIS depends on several Python packages. It is easiest to set up a separate conda environment to installed all required packages and to run COSMIS. All required packages can be installed when creating the conda environment, using the following command:
```bash
conda create --name cosmis --file requirements.txt
```
Obviously, you will need to install Miniconda or Anaconda before running the command above.

### Download required datasets
1. Get all transcript coding sequences from Ensembl.
```bash
wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz
```

2. Get the amino acid sequences of human reference proteome as annotated by UniProt.
```bash
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz
```

3. Download counts of unique variants at each amino acid position for all gnomAD annotated transcripts.

We preprocessed the gnomAD database and created a JSON formatted file that maps Ensembl stable transcript IDs to unique variant counts and variant types (missense or synonymous) for all position in the human proteome where a SNP variant was annotated by gnomAD. We have made this dataset available through [FigShare](https://figshare.com/ndownloader/files/31186919).

4. Get the mapping table from UniProt protein IDs to Ensembl stable transcript IDs.

A mapping from UniProt protein IDs to Ensembl stable transcript IDs is also required to run COSMIS. We have created such a mapping table and made it available through [FigShare](https://figshare.com/ndownloader/files/31186929)

5. Get transcript-level mutation probabilities.
   Transcript-level mutation probabilities are required to run COSMIS. You can get them from [FigShare](https://figshare.com/s/5f3e0fabc92a0ce59cdc). 
