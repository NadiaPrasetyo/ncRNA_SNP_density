# ncRNA_SNP_Density

## Overview
This repository contains a pipeline designed to analyze patterns of variation density in non-coding RNA (ncRNA). By utilizing various bioinformatics tools, this pipeline processes SNP (single nucleotide polymorphism) data to detect unusual patterns of variation in ncRNAs, with the ultimate goal of understanding how these variations might relate to ncRNA functions and their role in genomics.

The list of ncRNA sequences to be analyzed is provided in the `data/SNP-densities-and-RNA.csv` file.

## Directory Structure

- **.venv/**: Contains the virtual environment for Python dependencies. *This directory is ignored by Git and not included in the repository.*
- **bin/**: A collection of scripts for processing and analyzing data.
  - `bin/R/`: R scripts for statistical analysis and data visualization.
  - `bin/`: Additional Python, bash, and shell scripts for preprocessing, alignment, SNP processing, and result generation.
- **data/**: Input and output data for the analysis pipeline. *Note: Some large datasets (e.g., `data/Chromosomes_FASTA/`, `data/GEDS/`, `data/VCF/`) and certain files are excluded from the GitHub repository due to the `.gitignore` rules.*
  - `data/SNP-densities-and-RNA.csv`: Contains a list of ncRNAs to be analyzed.
  - Other subdirectories such as `data/Alignments/`, `data/VCF/`, etc., contain bioinformatics output data.
- **docs/**: Documentation and additional notes.
- **results/**: The output directory containing analysis results and visualizations.
- **.gitignore**: Specifies files and directories to be ignored by version control.
- **notes.txt**: Additional notes related to the project.

## Key Files in `bin/`
- `bin/R/Barplot_pangenome_variation.R`: R script to generate bar plots for variation across the pangenome.
- `bin/R/Enrichment_SNP155.R`: R script to analyze SNP enrichment in the pangenome.
- `bin/R/Z-score_normalise_pangenome_variation.R`: R script to normalize variation across the pangenome using Z-scores.
- `bin/align_geds.sh`: Shell script for aligning genome-edited data sets (GEDS).
- `bin/clean_fasta_CAPS.py`: Python script for cleaning FASTA sequences.
- `bin/collect_aligned.py`: Python script to collect aligned sequences.
- `bin/combine_fasta.py`: Python script to combine multiple FASTA files.
- `bin/fasta_to_fastq.py`: Converts FASTA format to FASTQ.
- `bin/fetch_fasta_NCBI.py`: Fetches FASTA sequences from NCBI.
- `bin/fetch_pangenome_data.py`: Fetches pangenome data from a specified source.
- `bin/fetch_snp_data.py`: Fetches SNP data for analysis.
- `bin/summarise_alignment.py`: Python script to summarize alignment data.
- `bin/summarise_pangenome_SNP.py`: Summarizes SNP data across the pangenome.
- `bin/summarise_SNP_data.py`: Summarizes SNP variation data.

## Workflow Overview
The pipeline follows several key steps to identify SNP variation density patterns in ncRNA:
1. **Data Collection**: Gather ncRNA sequences and associated SNP data from public databases and local files.
2. **Data Preprocessing**: Clean and filter FASTA sequences, align genomic data, and generate necessary files (e.g., VCF).
3. **Analysis**: Use statistical scripts (written in Python and R) to analyze SNP variation and enrichment in ncRNA sequences.
4. **Visualization**: Generate plots and summaries to visualize the distribution of SNPs across ncRNAs and identify patterns of variation density.
5. **Results**: Summarize findings in CSV files and produce visualizations (e.g., histograms, bar plots, UCSC links).

## Dependencies
- Python (3.x)
- R
- Required Python packages (installed via `requirements.txt`):
  - `bioservices 1.11.2`
  - `cattrs 24.1.2`
  - `certifi 2024.8.30`
  - `ipykernel 6.29.5`
  - `ipython 8.30.0`
  - `jupyter_client 8.6.3`
  - `jupyter_core 5.7.2`
  - `numpy 2.1.3`
  - `pandas 2.2.3`
  - `py_fasta_validator 0.6`
  - `pyBigWig 0.3.23`
  - `pybiomart 0.2.0`
  - `requests 2.32.3`
  - `seaborn 0.13.2`
  - `tqdm 4.67.1`
- R packages:
  - `BiocManager version 1.30.25`
  - `DescTools version 0.99.58`
  - `dplyr version 1.1.4`
  - `ggplot2 version 3.5.1`
  - `ggtext version 0.1.2`
  - `readr version 2.1.5`
  - `scales version 1.3.0`
  - `stringr version 1.5.1`
- External Libraries:
  - [GED-MAP](https://github.com/thomas-buechler-ulm/gedmap)
  - [SDSL library](https://github.com/simongog/sdsl-lite)


### Data
- Pangenome HPRC variant data: [UCSC download - HPRC decomposed.vcf.gz (899M)](https://hgdownload.soe.ucsc.edu/gbdb/hg38/hprc/)
- Human SNP variation data: [dbSnp155.bb (68G)](http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/)
- Human Chromosome sequences: [UCSC: hg38 chromosomes ](https://hgdownload.soe.ucsc.edu/goldenpath/hg38/chromosomes/)

### Installation
1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/ncRNA_SNP_Density.git
   cd ncRNA_SNP_Density

### Author
[Nadia Prasetyo](https://github.com/NadiaPrasetyo)
University of Otago, Department of Biochemistry
Research Associate

### License

This repository is licensed under the MIT License. See the LICENSE file for more information.
Acknowledgements

