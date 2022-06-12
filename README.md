# COSMIS
COSMIS is a novel framework for quantifying the 3D mutational constraint on amino acid sites in the human proteome. If you find COSMIS useful in your work, please consider citing the following paper: 
* Li, B., Roden, D.M., and Capra, J.A. [The 3D mutational constraint on amino acid sites in the human proteome](https://www.nature.com/articles/s41467-022-30936-x). *Nat. Commun.* **13**, 3273 (2022).

## What's in each folder?

### cosmis
The `cosmis` folder contains utility code that the top level cosmis code `cosmis.py`, `cosmis_batch.py`, and `cosmis_sp.py` depend on.

### figure-code-data
The `figure-code-data` folder contains standlone R code that can be run to reproduce figures in the main text and supplementary document of the manuscript.

### cosmis-scores
The `cosmis-scores` folder constains precomputed scores for all 16,533 proteins of the human reference proteome currently covered by the framework.

### scripts
The `scripts` folder contains all scripts written in Python that were called to obtain processed datasets in the `data` folder.

### structures
The `structures` folder contains all protein structure in the PDB format based on which COSMIS scores were computed.

### supplementary-data
The `supplementray-data` folder contains all supplementary tables referred to in the published COSMIS paper.

## Using the COSMIS framework
It is recommended that interested users of the COSMIS framework download precomputed COSMIS scores from this repository. However, should you need to run COSMIS using custom-built protein structural models, or to compute COSMIS scores based on protein-protein complexes, please follow the following steps.

### Clone COSMIS
Clone COSMIS to a local directory.
```bash
git clone https://github.com/CapraLab/cosmis.git
```
Note that cloning might fail as this repository is tracked with [Git Large File Storage](https://git-lfs.github.com/) and is over the data quota currently allowed by Git LFS. Please check later as we sorting out this quota issue.

### Set up conda environment
COSMIS depends on several Python packages. It is easiest to set up a separate conda environment to installed all required packages and to run COSMIS. All required packages can be installed when creating the conda environment, using the following commands:
```bash
# if your platform is Linux based, run this
conda create --name cosmis --file requirements_linux.txt

# if your platform is OSX based, run this
conda create --name cosmis --file requirements_osx.txt

# activate the environment
conda activate cosmis

# then also use pip to install wget under the environment
pip install wget
```
Obviously, you will need to install Miniconda or Anaconda before running the command above.

### Download required datasets
Some of the required datasets are already made available within this repository (in the `database_files` folder). However, due to limits on file size, we had to made larger files available through other means. Follow the following steps to get all required input datasets.
1. Get all transcript coding sequences from Ensembl (this dataset is already available in `database_files`).
```bash
wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz
```

2. Get the amino acid sequences of human reference proteome as annotated by UniProt (this dataset is already available in `database_files`).
```bash
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz
```

3. Download counts of unique variants at each amino acid position for all gnomAD annotated transcripts.

   We preprocessed the gnomAD database and created a JSON formatted file that maps Ensembl stable transcript IDs to unique variant counts and variant types (missense or synonymous) for all position in the human proteome where a SNP variant was annotated by gnomAD. We have made this dataset available through [FigShare](https://figshare.com/ndownloader/files/31186919).

4. Get the mapping table from UniProt protein IDs to Ensembl stable transcript IDs (this dataset is already available in `database_files`).

   A mapping from UniProt protein IDs to Ensembl stable transcript IDs is also required to run COSMIS. We have created such a mapping table and made it available through [FigShare](https://figshare.com/ndownloader/files/31186929)

5. Get transcript-level mutation probabilities (this dataset is already available in `database_files`).
   
   Transcript-level mutation probabilities are required to run COSMIS. You can get them from [FigShare](https://figshare.com/s/5f3e0fabc92a0ce59cdc). 

### Run COSMIS
Depending on whether you want run COSMIS on a single monomeric protein or homo-multimeric protein, or a list of monomeric proteins, the script and setup are slightly different.

1. Run COSMIS on a single monomeric or homo-multimeric protein.

1.1 Use the following JSON formatted template to supply paths to database files.
```bash
{
    "ensembl_cds": "/path/to/Homo_sapiens.GRCh38.cds.all.fa.gz",
    "uniprot_pep": "/path/to/UP000005640_9606.fasta.gz",
    "gnomad_variants": "/path/to/gnomad_filtered/gnomad_variant_counts_hg38.json",
    "uniprot_to_enst": "/path/to/uniprot_to_enst.json",
    "enst_mp_counts": "/path/to/mutation_probs.tsv"
}
```
Then, save the file as `data_paths.json`, for example.

1.2 If you'd like to compute COSMIS score WITHOUT accounting for contacts from neighboring subunits. Run this command
```bash
python cosmis_sp.py -c data_paths.json -u <UniProt ID> -p <PDB file> --chain <chain ID of subunit> -o monomeric_cosmis.tsv
```

1.3 If you'd to compute COSMIS score accounting for contacts from neighboring subunits. Add `--multimer` to the command above, i.e.
```bash
python cosmis_sp.py -c data_paths.json -u <UniProt ID> -p <PDB file> --chain <chain ID of subunit> -o multimeric_cosmis.tsv --multimer
```

1.4 One can also run the following command to compute COSMIS scores for a subunit which is part of a hetero-oligomeric protein complex. However, `cosmis_complex.py` has not been thoroughly tested or benchmarked. Interpret the results with caution and let us know if you find anything buggy.
```bash
python cosmis_complex.py -c data_path.json -i <chain_to_uniprot_mapping file> -p <PDB file> --chain <chain ID of subunit> -o <output file>

# example
cd examples/
python ../cosmis_complex.py -c cosmis_config.json -i KCNQ1_chain_to_uniprot.txt -p KCNQ1.pdb -o KCNQ1_cosmis.tsv --chain A
```

2. Run COSMIS on a list of monomeric proteins whose structures were obtained from AlphaFold database or SWISS-MODEL repository.

2.1 Use the following JSON formatted template to supply paths to database files.
```bash
{
    "ensembl_cds": "/path/to/Homo_sapiens.GRCh38.cds.all.fa.gz",
    "uniprot_pep": "/path/to/UP000005640_9606.fasta.gz",
    "gnomad_variants": "/path/to/gnomad_filtered/gnomad_variant_counts_hg38.json",
    "uniprot_to_enst": "/path/to/uniprot_to_enst.json",
    "enst_mp_counts": "/path/to/mutation_probs.tsv"
    "pdb_dir": "/path/to/pdb_files/",
    "output_dir": "./"
}
```
Then, save the file as `data_paths.json`, for example.

2.2 If the structures were obtained from AlphaFold database, run the following command
```bash
python cosmis_batch.py -c data_paths.json -i <input.txt> -d AlphaFold -l af_cosmis.log
```
The input file `input.txt` contains on each line a pair of UniProt ID and PDB filename. For example
```bash
A0A024R1R8 A0/A0/AF-A0A024R1R8-F1-model_v1.pdb
A0A024RBG1 A0/A0/AF-A0A024RBG1-F1-model_v1.pdb
A0A024RCN7 A0/A0/AF-A0A024RCN7-F1-model_v1.pdb
A0A075B6H5 A0/A0/AF-A0A075B6H5-F1-model_v1.pdb
A0A075B6H7 A0/A0/AF-A0A075B6H7-F1-model_v1.pdb
A0A075B6H8 A0/A0/AF-A0A075B6H8-F1-model_v1.pdb
A0A075B6H9 A0/A0/AF-A0A075B6H9-F1-model_v1.pdb
A0A075B6I0 A0/A0/AF-A0A075B6I0-F1-model_v1.pdb
A0A075B6I1 A0/A0/AF-A0A075B6I1-F1-model_v1.pdb
A0A075B6I3 A0/A0/AF-A0A075B6I3-F1-model_v1.pdb
```
assuming that the base directory where the PDB files are stored is `/path/to/pdb_files/`.

2.3 If the structures were obtained from SWISS-MODEL repository, run the following command 
```bash
python cosmis_batch.py -c data_paths.json -i <input.txt> -d SWISS-MODEL -l swiss_model_cosmis.log
```
The input file `input.txt` contains on each line a pair of UniProt ID and PDB filename. For example
```bash
A8MWA4 A8/MW/A4/swissmodel/109_299_5v3m.1.C.pdb
A8MWD9 A8/MW/D9/swissmodel/3_76_4wzj.3.G.pdb
A8MWL7 A8/MW/L7/swissmodel/11_110_2loo.1.A.pdb
A8MX76 A8/MX/76/swissmodel/21_683_3bow.1.A.pdb
A8MXE2 A8/MX/E2/swissmodel/68_331_7jhi.1.A.pdb
A8MXQ7 A8/MX/Q7/swissmodel/136_377_4qfv.1.A.pdb
A8MXT2 A8/MX/T2/swissmodel/95_311_2wa0.1.A.pdb
A8MXU0 A8/MX/U0/swissmodel/23_60_2lwl.1.A.pdb
A8MXY4 A8/MX/Y4/swissmodel/200_532_5v3m.1.C.pdb
A8MYX2 A8/MY/X2/swissmodel/110_172_5cwg.1.A.pdb
```
again, assuming that the base directory where the PDB files are stored is `/path/to/pdb_files/`.
