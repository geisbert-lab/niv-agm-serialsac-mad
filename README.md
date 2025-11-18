# Transcriptomic analysis for *Progression of Nipah virus infection in African green monkeys: role of macrophages and dendritic cells in early virus spread*

> This manuscript was submitted for review to Cell Reports (`CELL-REPORTS-D-25-07505`) on 17-Nov-25.

## Data availability

Raw files are available via NCBI GEO.
- **[`GSE309971`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309971):** Nanostring nCounter RCC files 
- **[`GEO TBD`](linktbd):** RNA-seq FASTQ files

Count matrices are available via NCBI GEO or in this repo's `data/` directory.

Code to generate figures is in `analysis.r`, and output from this script are in `analysis/`.

## Methods and preprocessing

### Targeted transcriptomics of circulating blood (Nanostring)

The expression of 770 host mRNAs was quantified from whole blood RNA via the [Nanostring NHP Immunology v2 panel](https://nanostring.com/products/ncounter-assays-panels/ncounter-rna-assays/immunology/nhp-immunology/) ([PMID 29116224](https://pubmed.ncbi.nlm.nih.gov/29116224/); Bruker cat. #115000276) according to the manufacturer’s instructions. Raw RCC files were loaded into [nSolver v4.0](https://nanostring.com/products/ncounter-analysis-system/ncounter-analysis-solutions/), and background thresholding was performed using the default parameters. Samples that failed the nSolver internal quality checks were removed. Thresholded count matrices were exported from nSolver and analyzed with [limma v3.62.1](https://bioconductor.org/packages/release/bioc/html/limma.html) ([edgeR v4.4.1](https://bioconductor.org/packages/release/bioc/html/edgeR.html)) in R v4.5.0. Samples were binned by days postinfection (3, 4, and 5 DPI) and compared to baseline (0 DPI). mRNAs with an FDR-adjusted p-value and a log2 fold change >1 or <-1 were considered significantly differentially expressed. Canonical signaling pathways, functions, and upstream regulators were identified via [Ingenuity Pathway Analysis](https://digitalinsights.qiagen.com/products-overview/discovery-insights-portfolio/analysis-and-visualization/qiagen-ipa) ([PMID 24336805](https://pubmed.ncbi.nlm.nih.gov/24336805/)). 

## RNA-seq of tonsil and lung

RNA from tonsil and lung samples homogenized in RLT was used to construct poly-A-selected sequencing libraries. Samples were sequenced to a target depth of 50 million 150x150 paired end reads on the NextSeq 550 (Illumina). Raw FASTQ files were aligned to the combined AGM ([ChlSab1, Ensembl v112](https://useast.ensembl.org/Chlorocebus_sabaeus/Info/Index)) and NiV Bangladesh ([GCA_003147825.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_003147825.1/)) genomes using STAR v2.7.11b ([PMID 23104886](https://pubmed.ncbi.nlm.nih.gov/23104886/)). PCR duplicates were removed with samtools v1.20 ([PMID 19505943](https://pubmed.ncbi.nlm.nih.gov/19505943/)), and count matrices were generated with featureCounts ([PMID 24227677](https://pubmed.ncbi.nlm.nih.gov/24227677/)) via Rsubread v2.20.0 ([PMID 30783653](https://pubmed.ncbi.nlm.nih.gov/30783653/)) in R v4.1.3. 


These processing steps were completed via the [Stampede3 supercomputer](https://tacc.utexas.edu/systems/stampede3/). The STAR index was created with `index-star.sh`. Code for alignment and PCR deduplication can be found in `alignment-star.sh`, and quantification code is in `quantification-featurecounts.sh`. Output logs for each job can be found in `logs/`.

Differential expression analyses were performed with DESeq2 v1.49.1 ([PMID 25516281](https://pubmed.ncbi.nlm.nih.gov/25516281/)) in R v4.5.0 by comparing each tissue and time point postinfection to baseline. Genes with an FDR-adjusted p-value and a log2 fold change >2 or <-2 were considered significantly differentially expressed and used as input for [Ingenuity Pathway Analysis](https://digitalinsights.qiagen.com/products-overview/discovery-insights-portfolio/analysis-and-visualization/qiagen-ipa) ([PMID 24336805](https://pubmed.ncbi.nlm.nih.gov/24336805/)). The code for these analyses and figure generation can be found in `analysis.r`.
