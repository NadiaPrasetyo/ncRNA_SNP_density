# ncRNA_SNP_Density

## Overview
This repository contains a pipeline designed to analyze patterns of variation density in non-coding RNA (ncRNA). 
By utilizing various bioinformatics tools, this pipeline processes SNP (single nucleotide polymorphism), small indels,
and other variation data to detect unusual patterns of variation in ncRNAs, with the ultimate goal of understanding
how these variations might relate to ncRNA functions and their role in genomics.

The list of ncRNA sequences to be analyzed is provided in the `data/SNP-densities-and-RNA.csv` file.

## Directory Structure

- **bin/**: A collection of scripts for processing and analyzing data.
  - `bin/R/`: R scripts for statistical analysis and data visualization.
  - `bin/`: Additional Python and shell scripts for preprocessing, alignment, SNP processing, and result generation.
- **data/**: Input and output data for the analysis pipeline. *Note: Some large datasets must be downloaded external to the github.*
  - `data/SNP-densities-and-RNA.csv`: Contains a list of ncRNAs to be analyzed and their SNP density as detected by SNP finding algorithms.
  - Other subdirectories such as `data/Alignments/`, `data/VCF/`, etc., contain fetch output data that is too large for tthe github repository.
  - `data/FASTA`: sequences of each of the ncRNA as taken from NCBI.
  - `data/fastq`: fastQ converted sequence of each ncRNA for analysis with GED-MAP.
- **docs/**: Documentation and additional notes.
  -  `docs/Map-notes.ipynb`: Comparison between GED-MAP alignment and BLASTn hit alignments of each ncRNA against the human chromosome reference/construct.
- **results/**: The output directory containing analysis results and visualizations.
  - `results/UCSC_tracks_100x_ZoomOut`: UCSC tracks downloaded for each ncRNA to reveal pangenome variation in a broader area.
- **.gitignore**: Specifies files and directories to be ignored by version control.
- **notes.txt**: Additional notes related to the project.

## Key Files in `bin/`

- `bin/R/Barplot_pangenome_variation.R`: R script to generate bar plots for variation across the pangenome.
- `bin/R/Enrichment_SNP155.R`: R script to analyze SNP enrichment in the pangenome.
- `bin/R/Z-score_normalise_pangenome_variation.R`: R script to normalize variation across the pangenome using Z-scores.
- `bin/align_each_genome.sh`: Shell script for aligning individual genomes.
- `bin/align_geds.sh`: Shell script for indexing and aligning pangenome construct data sets (GEDS).
- `bin/clean_fasta_CAPS.py`: Python script for cleaning FASTA sequences, capitalizing all of the bases (ACGTN).
- `bin/collect_aligned.py`: Python script to collect aligned sequences.
- `bin/combine_fasta.py`: Python script to combine multiple FASTA files.
- `bin/count_num_var.py`: Python script to count the number of variations in the data.
- `bin/EXTENDED_fetch_pangenome_data.py`: Python script to fetch extended pangenome data.
- `bin/EXTENDED_fetch_snp_data.py`: Python script to fetch extended SNP data.
- `bin/fasta_to_fastq.py`: Python script to convert FASTA format to FASTQ.
- `bin/fetch_animals_loci.py`: Python script to fetch loci information for different animals.
- `bin/fetch_fasta_NCBI.py`: Python script to fetch FASTA sequences from NCBI.
- `bin/fetch_fasta_seq.py`: Python script to fetch FASTA sequences.
- `bin/fetch_pangenome_data.py`: Python script to fetch pangenome data from UCSC decomposed.vcf.
- `bin/fetch_snp_data.py`: Python script to fetch SNP data from UCSC dbSnp155.bb.
- `bin/filter_vcf_by_genome.py`: Python script to filter VCF files by genome.
- `bin/gene_range_calc.py`: Python script to calculate gene ranges for analysis.
- `bin/generate_UCSC_link.py`: Python script to generate UCSC genome browser links.
- `bin/parse_chrom_loc.py`: Python script to parse chromosome locations.
- `bin/parse_each_genome.sh`: Shell script for parsing individual genomes.
- `bin/process_blast_result.py`: Python script to process BLAST results.
- `bin/process_near_perf_maps.py`: Python script to process near-perfect maps.
- `bin/rfam_cm_animals.py`: Python script to fetch Rfam CM data for animals.
- `bin/run_gedmap.sh`: Shell script to run GEDMAP analysis.
- `bin/separate_genome_alignment.py`: Python script to separate genome alignment data.
- `bin/separate_vcf.py`: Python script to separate VCF data.
- `bin/sort_and_clean_vcf.py`: Python script to sort and clean VCF files.
- `bin/SPKM_count.py`: Python script to count SPKM (Scaled Pseudo-Kmers).
- `bin/summarise_alignment.py`: Python script to summarize alignment data.
- `bin/summarise_pangenome_SNP.py`: Python script to summarize SNP data across the pangenome to plottable data.
- `bin/summarise_SNP_data.py`: Python script to summarize SNP variation data to plottable data.

## Workflow Overview
The pipeline follows several key steps to identify SNP variation density patterns in ncRNA:
1. **Data Collection**: Gather ncRNA sequences and associated SNP data from public databases and local files.
2. **Data Preprocessing**: Clean and filter FASTA sequences, align genomic data, and generate necessary files (e.g., VCF).
3. **Analysis**: Use statistical scripts (written in Python and R) to analyze SNP variation and enrichment in ncRNA sequences.
4. **Visualization**: Generate plots and summaries to visualize the distribution of SNPs across ncRNAs and identify patterns of variation density.
5. **Results**: Summarize findings in CSV files and produce visualizations (e.g., histograms, bar plots, UCSC links).

## GED-MAP analysis of a pangenome construct of ncRNA to human chromosome sequence.
1. Download the VCF data of pangenome variation in the body - UCSC decomposed.vcf.
2. Separate the VCF to different chromosomes to be applied individually. Done by `bin/separate_vcf.py`. Then sort and clean the variation according to `bin/sort_and_clean_vcf.py`.
3. Download chromosome sequence data from NCBI (RefSeq) - preferably at least chr1 - chr22.
4. Use GED-MAP (installed separately) to parse the chromosome sequence and variation for the respective chromosome to create a pangenome chromosome construct. Done in `bin/run_gedmap.sh`.
5. Use GED-MAP to index and align each ncRNA fastq sequences with the chromosome pangenome construct to find any alignments and hits where the ncRNA sequence is found in the chromosome. Done in `bin/align_geds.sh`.
6. Collect the aligned results (SAM file) using `bin/collect_aligned.py`.
7. Summarise the aligned results using `bin/summarise_alignment.py`.And analyse in the docs `docs/Map-notes.ipynb`.

## Key Files in `data/`

- `data/Alignments/`: Contains alignment data for various genomic sequences. This directory is not tracked due to ists large size but the contents can be made using the scripts in bin/.
- `data/animals/`: Directory containing animal-related genomic data. This directory is not tracked due to ists large size but the data can be downloaded separately
  - `data/animals/output/`: Directory for output data of animal-rfam cmscan. Contains the deoverlapped tblout scan of each animal genome with the different rfam models.
  - `data/animals/rfam/`: Directory for Rfam models for animals.
  - `data/animals/bosTau9.fa`: FASTA file for Bos taurus (cow) genome.
  - `data/animals/bosTau9.vcf`: VCF file for Bos taurus (cow) variants.
  - `data/animals/danRer11.fa`: FASTA file for Danio rerio (zebrafish) genome.
  - `data/animals/danRer11.vcf`: VCF file for Danio rerio (zebrafish) variants.
  - `data/animals/galGal6.fa`: FASTA file for Gallus gallus (chicken) genome.
  - `data/animals/galGal6.vcf`: VCF file for Gallus gallus (chicken) variants.
  - `data/animals/hg38.fa`: FASTA file for Homo sapiens (human) genome.
  - `data/animals/hg38.vcf`: VCF file for Homo sapiens (human) variants.
  - `data/animals/mm39.fa`: FASTA file for Mus musculus (mouse) genome.
  - `data/animals/mm39.vcf`: VCF file for Mus musculus (mouse) variants.
  - `data/animals/rn7.fa`: FASTA file for Rattus norvegicus (rat) genome.
  - `data/animals/rn7.vcf`: VCF file for Rattus norvegicus (rat) variants.
- `data/BLASTN_RefSeq/`: Directory containing BLASTN results from NCBI RefSeq downloaded as separate JSONs. This directory is not tracked due to its large size.
- `data/FASTA/`: Directory containing FASTA sequence files for various ncRNAs.
- `data/fastq/`: Directory containing FastQ converted sequences for analysis with GED-MAP.
- `data/GEDS/`: Directory containing data for genomic datasets related to the pangenome. This directory is not tracked due to its large size but the contents can be made using the scripts in bin/.
- `data/separated_genome_sam/`: Directory containing separated SAM files for genome data.
- `data/VCF/`: Directory containing genome VCF files. Untracked due to size, contents can be made using scripts in bin/.
- `data/adjusted_genome_alignment.sam`: SAM file for adjusted genome alignment.
- `data/aligned_genome_reads.sam`: SAM file for aligned genome reads.
- `data/aligned_reads.sam`: SAM file for aligned sequence reads.
- `data/animals_var.csv`: CSV file containing variant data for different animals.
- `data/BLASTn_results_summary.csv`: CSV file summarizing BLASTn results.
- `data/combined_sequences.fa`: Combined FASTA file for sequences across multiple sources.
- `data/EXTENDED_filtered_variants.csv`: CSV file containing extended filtered variant data.
- `data/EXTENDED_pangenome_summary.csv`: CSV file summarizing extended pangenome data.
- `data/EXTENDED_SNP_data.txt`: Text file containing extended SNP data.
- `data/filtered_ncRNA_list.csv`: CSV file containing a filtered list of ncRNAs.
- `data/filtered_variants.csv`: CSV file containing filtered variants.
- `data/filtered_variants.vcf`: VCF file containing filtered variants.
- `data/other_animals_RNA.csv`: CSV file containing RNA data for other animals.
- `data/pangenome_summary.csv`: CSV file summarizing pangenome data.
- `data/SNP_data.txt`: Text file containing SNP data for the pangenome.
- `data/SNP-densities-and-RNA.csv`: CSV file containing SNP density data and corresponding RNA sequences.

## Dependencies
- Python (3.10.12)
- R (4.4.2)
- Required Python packages:
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
  - [Rfam CM](https://docs.rfam.org/en/latest/genome-annotation.html)
  - [UCSC twoBitToFa (9.0M)](https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/)


### Data
- Pangenome HPRC variant data: [UCSC download - HPRC decomposed.vcf.gz (899M)](https://hgdownload.soe.ucsc.edu/gbdb/hg38/hprc/)
- Human SNP variation data: [dbSnp155.bb (68G)](http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/)
- Human Chromosome sequences: [UCSC: hg38 chromosomes ](https://hgdownload.soe.ucsc.edu/goldenpath/hg38/chromosomes/)
- Human genome: [hg38.2bit (797M)](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)
- Chicken genome: [galGal6.2bit (299M)](https://hgdownload.soe.ucsc.edu/goldenPath/galGal6/bigZips/)
- Chicken variation: [9031_GCA_000002315.5_current_ids.vcf.gz (361M)](https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_6/by_species/gallus_gallus/GCA_000002315.5/)
- Mouse genome: [mm39.2bit (681M)](https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/)
- Mouse variation: [10090_GCA_000001635.9_current_ids.vcf.gz (1.6G)](https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_6/by_species/mus_musculus/GCA_000001635.9/)
- Rat genome: [rn7.2bit (660M)](https://hgdownload.soe.ucsc.edu/goldenPath/rn7/bigZips/)
- Rat variation: [10116_GCA_015227675.2_current_ids.vcf.gz (100M)](https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_6/by_species/rattus_norvegicus/GCA_015227675.2/)
- Zebrafish genome: [danRer11.2bit (421M)](https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/)
- Zebrafish variation: [7955_GCA_000002035.4_current_ids.vcf.gz (201M)](https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_6/by_species/danio_rerio/GCA_000002035.4/)
- Cow genome: [bosTau9.2bit(680M)](https://hgdownload.soe.ucsc.edu/goldenPath/bosTau9/bigZips/)
- Cow variation: [9913_GCA_002263795.2_current_ids.vcf.gz (1.5G)](https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_6/by_species/bos_taurus/GCA_002263795.2/)

### Author
[Nadia Prasetyo](https://github.com/NadiaPrasetyo)
University of Otago, Department of Biochemistry
Research Associate

### License

This repository is licensed under the MIT License. See the LICENSE file for more information.
Acknowledgements

