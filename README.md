# RNA‑seq differential expression pipeline for <i> Anopheles stephensi </i> larvae

## A small, end‑to‑end toolkit to quantify transcripts with Salmon, build a custom org.*.eg.db from a GAF file, and run differential expression plus GO/KEGG enrichment in R.

### It includes:

PBS/Apptainer scripts to index a transcriptome and quantify paired‑end reads with Salmon.

A Python utility to parse FASTA headers and produce a tx2gene/annotation CSV.

An R script to build a custom OrgDb package (org.Astephensi.eg.db) from a GAF 2.2 file.

An R workflow for DESeq2, VST/PCA, KEGG ORA, GO GSEA, and GO ORA.


The R script expects UCI larval stages coded as IS1/IS2/IS3/IS4 and a quant.sf per sample.

# Requirements
## System

Linux HPC (tested with PBS/Torque) and Apptainer for containerized Salmon runs.

Disk space for transcriptome indexing and quant outputs.

## Software

Salmon (inside your Apptainer image).

Python 3.8+ (standard library only) for tx2gene_creator.py.

R ≥ 4.5 with packages:

DESeq2, tximport, apeglm, pheatmap

dplyr, tibble, ggplot2, readr

clusterProfiler, enrichplot

AnnotationForge, AnnotationDbi, biomaRt, GO.db

(after build) org.Astephensi.eg.db

# Quick start
0) Set up your folders

Create a project directory with subfolders like:

RNAseq_EthiopiaNUCI/
  Ref_transcriptome/
  Fastq_files/UCI_lab/
  quant_files/UCI_quants/


## Prepare:

Anstep_UCI_V1.0_rna.fna(.gz) in Ref_transcriptome/

UCI_ANSTEP_V1.0_gene_ontology.gaf in Ref_transcriptome/

Paired FASTQs in Fastq_files/UCI_lab/

A file UCIlab_samples.txt listing sample IDs (without _1/_2 suffix)

## 1) Index the transcriptome with Salmon (HPC)

Edit paths in hpc/indexing_reference_salmon.sh to match your environment and submit to PBS:

qsub hpc/indexing_reference_salmon.sh


This script runs salmon index against Anstep_UCI_V1.0_rna.fna(.gz) and writes Anstep_UCI_V1.0_rna_index.

## 2) Quantify reads with Salmon (HPC)

Edit paths in hpc/quant_salmon.sh, ensure UCIlab_samples.txt is present, then:

qsub hpc/quant_salmon.sh


The script loops over sample IDs and runs salmon quant (library autodetect, selective alignment, GC‑bias correction) producing one quant.sf per sample under quant_files/UCI_quants/.

## 3) Create transcript annotations (tx2gene source)

From Reference_transcriptome/:

python3 python/tx2gene_creator.py


This parses FASTA headers and writes Anstep_UCI_V1.0_transcripts_info.csv with:
accession, predicted, organism, product, gene_id (e.g., LOC...), molecule.

The DE script later trims transcript version suffixes (e.g., .1) and removes the LOC prefix from gene IDs for KEGG/GO compatibility.

## 4) Build the custom OrgDb from GAF (optional but recommended)

From Reference_transcriptome/:

source("R/gaf2orgDB.R")


This will:

Read a GAF 2.2 file, drop NOT qualifiers by default,

Construct gene_info and GO tables,

Build & install org.Astephensi.eg.db (version set in the script),

Validate with keytypes() / columns() checks.

Make sure maintainer and author strings follow the script’s required format.

## 5) Differential expression & enrichment analysis in R

Update paths at the top of R/Diff_Exp_UCI.R, then run:

source("R/Diff_Exp_UCI.R")


# What it does (high level):

Input: UCI_sample_info.csv (columns: sampleName, stage, quant_dir), builds file list of quant.sf.

Import: tximport aggregates Salmon transcript counts to genes using tx2gene derived from Anstep_UCI_V1.0_transcripts_info.csv.

DESeq2: Builds four DESeqDataSets (ref levels IS1...IS4), prefilters low counts, fits models, and shrinks LFC with apeglm.

Comparisons: Extracts IS2 vs IS1, IS3 vs IS2, IS4 vs IS3 with the correct reference models (pay attention to reference level when interpreting up/down sets).

QC/Exploration: VST → PCA plot; sample distance heatmap.

Exports: Significant up/down gene tables (padj < 0.05 & |log2FC| > 1).

KEGG ORA: Uses KEGG organism code aste for An. stephensi; saves barplots & dotplots.

GO GSEA: Uses the custom org.Astephensi.eg.db to assemble TERM2GENE and TERM2NAME (BP only), runs clusterProfiler::GSEA, and saves dotplots, ridgeplots.

GO ORA: Enrichment of significant up/down sets against the tested gene universe, with dotplots and barplots.

Outputs are saved under ./UCI_Diff_Exp_plots/ and (optionally) CSVs for normalized counts and DE results.

# Data conventions & tips

Stages: Ensure your stage factor is ordered IS1, IS2, IS3, IS4 in the sample table.

Paths: The scripts are path‑sensitive; search for setwd(...), cd, and hard‑coded directories and update as needed.

KEGG IDs: Script removes LOC from gene IDs to match KEGG’s ncbi-geneid style.

GAF 2.2: If your GAF schema differs, adjust column handling in gaf2orgDB.R.

# Troubleshooting

Reference levels & contrasts: When asking “genes up in A vs B,” the model must use B as reference to interpret signs correctly. Double‑check your relevel() calls and results(..., contrast=...).

Missing KEGG/GO hits: Confirm ID types (Entrez/NCBI‑GeneID) and organism code (aste) and that LOC prefixes are removed where required.

PBS/Apptainer: If apptainer is not on $PATH on compute nodes, load the appropriate module or use an absolute path; verify the SIF image path in both HPC scripts.

# Citations / Acknowledgments

Patro et al., Salmon: fast transcript quantification (Nat Methods 2017).

Love et al., DESeq2: differential analysis of count data (Genome Biol 2014).

Yu et al., clusterProfiler: enrichment analysis (OMICS 2012; updates).

AnnotationForge / AnnotationDbi / GO.db (Bioconductor).

Apptainer for containerized HPC runs.



# Maintainer

Isuru Gunarathna — Baylor University
For questions or contributions, please open an issue or pull request.
