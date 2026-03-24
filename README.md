# CORAL
RShiny app built for RNASeq DEG result filtering, heatmaps, and GO term visualization 

Run this app locally in RStudio to use

## Setup
Before running, please ensure the necessary packages are installed using `CORAL_requirements.Rmd`

## Input
Three inputs are required (all files must be either csv or tsv format):

-A gene count (TPM) table, or other matrix of abundance data per sample, where each row is a unique gene and each column is an individual sample

-A metadata table, where the first column contains sample IDs matching the column names of the gene counts, and the following columns are independent treatment variables

-EITHER a folder containing DESeq2 output files (they must contain gene names matching those in the gene count table, and have at least an adjusted p-value and log2foldchange column) to be filtered for significant DEGs

-OR a single file containing a list of unique genes (this gene list cannot be filtered, please manually note the padj and log2fc values)

## Usage
To save any file, you must select a file in the `Data processing` tab. Here you can also choose the padj and log2fc thresholds

Heatmaps will be generated using the currently selected padj and log2fc values. Modifying them will actively change existing figures

## GO term plots
To generate GO term plots, please input two files (they must be in csv or tsv format, txt files will be accepted)

### GO <-> Gene
This file must contain 2 columns (in order):
```
GO ID number    Gene name (as found in gene table)
```

### GO <-> Term
This file must contain 2 columns (in order):
```
GO ID number    GO term description
```

GO analyses will be plotted based on the corresponding GO IDs and related terms