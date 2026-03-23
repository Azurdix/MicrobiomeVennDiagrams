# Core Microbiome Venn Diagrams - PRD / PTSD / Compartments

## Overview

This script identifies and visualizes the **core microbiome** of two spider species:

- **PRD** – *Pardosa lugubris*
- **PTSD** – *Parasteatoda tepidariorum*

The analysis compares bacterial taxa detected in three sample compartments:

- **ENV** – environmental samples
- **SILK** – silk-associated samples
- **EGGS** – egg-associated samples

Core taxa are identified separately for each species-compartment combination and then compared using **Venn diagrams**. The script also exports supplementary tables listing shared and unique taxa for each comparison.

---

## Purpose

The goal of this script is to detect bacterial families or genera that are:

- unique to a specific compartment,
- shared between selected compartments,
- shared across the full **ENV → SILK → EGGS** gradient,
- comparable between the two spider species.

This provides a simple and publication-friendly way to explore overlap patterns within the **core microbiome**.

---

## Input data

The script requires three tab-separated input files:

### 1. Metadata file
Path:
`01_Metadata/metadata_plik.tsv`

Required columns:

- `sample-id`
- `group` — species label (`PRD` / `PTSD`)
- `type` — compartment label (`ENV` / `SILK` / `EGGS`)

### 2. Feature table
Path:
`06_Exports/deseq2_input/feature-table.tsv`

This should be a QIIME2-exported feature table containing ASV counts per sample (mitochondria + chloroplasts filtered out).

### 3. Taxonomy file
Path:
`06_Exports/deseq2_input/taxonomy_export/taxonomy.tsv`

This file should contain taxonomy assignments for each feature (ASV), exported from QIIME2.

---

## How the script works

The workflow consists of the following main steps:

1. **Load input data**
   - metadata
   - feature table
   - taxonomy table

2. **Extract taxonomy**
   - family names are extracted from taxonomy strings
   - genus names are extracted from taxonomy strings

3. **Match valid samples**
   - only samples present in both metadata and feature table are retained

4. **Aggregate ASV counts**
   - counts are summed to the selected taxonomic level:
     - `family`
     - `genus`

5. **Define core microbiome**
   - a taxon is considered part of the core microbiome if it:
     - has more than `10` reads in a sample,
     - is present in at least `66%` of samples within a given group

6. **Filter taxonomy labels**
   - ambiguous or non-informative labels are removed, such as:
     - `Unassigned`
     - `Unknown`
     - `uncultured`
     - unresolved subgroup-like labels
     - overly broad taxonomic placeholders

7. **Generate Venn diagrams**
   - two-set and three-set comparisons are produced

8. **Save outputs**
   - plots in multiple formats
   - CSV tables with shared and unique taxa

---

## Core microbiome definition

The default thresholds used in this script are:

- **Count threshold:** `> 10` reads per sample
- **Prevalence threshold:** `>= 66%` of samples within a group

These values can be modified in the **Settings** section of the script.

---

## Output files

The script creates the following output directories:

- `plots/` — graphical outputs
- `supplementary_tables/` — summary and overlap tables
- `raw_sets/` — exported raw taxon sets

### Plot formats
Each Venn diagram is saved as:

- `.png`
- `.pdf`
- `.svg`

### Supplementary tables
For each comparison, the script exports CSV files containing:

- taxa shared between all sets
- taxa shared between selected pairs only
- taxa unique to each set
- summary counts for all Venn regions

---

## Main settings

The most important user-defined parameters are:

### Taxonomic level
```r
TAX_LEVEL <- "family"   # "genus" or "family"
