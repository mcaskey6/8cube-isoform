# 8cube isoform analysis

A Snakemake workflow for single-nucleus isoform analysis of the IGVF **8-cube** mouse dataset. It combines Oxford Nanopore long reads and Illumina short reads from the same Parse Biosciences (SPLiT-seq) libraries to build an isoform-level transcriptome on a C57BL/6 telomere-to-telomere (T2T) genome, then quantifies those isoforms in each cell.

## Background and data

The 8-cube dataset profiles eight tissues from each of the eight Collaborative Cross founder strains with single-nucleus RNA-seq, as part of the [IGVF consortium](https://igvf.org/):

> Rebboah E, Weber R, *et al.* Systematic cell-type resolved transcriptomes of 8 tissues in 8 lab and wild-derived mouse strains capture global and local expression variation. *Cell Genomics* 6(4) (2026). <https://doi.org/10.1016/j.xgen.2025.101108>

- **Short reads (Illumina):** public on the [IGVF Data Portal](https://data.igvf.org/). The paper's Table S2 links the raw FASTQ measurement set for every experiment and subpool. For example, experiment `igvf_010` (kidney, gastrocnemius) subpool 2 is [IGVFDS4866EOBI](https://data.igvf.org/measurement-sets/IGVFDS4866EOBI/).
- **Long reads (Nanopore):** not yet published. A subset was used in the lr-kallisto paper:

  > Loving RK, Sullivan DK, Reese F, Rebboah E, *et al.* Long-read sequencing transcriptome quantification with lr-kallisto. *PLoS Computational Biology* (2025). <https://doi.org/10.1371/journal.pcbi.1013692>

The barcode extraction steps follow the lr-kallisto example notebook in [pachterlab/LSRRSRLFKOTWMWMP_2024](https://github.com/pachterlab/LSRRSRLFKOTWMWMP_2024).

## What the workflow does

1. **Reference:** downloads the C57BL/6 T2T genome ([GCA_056825265.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_056825265.1/), chr1–19, X, Y, MT) and lifts the GENCODE vM39 (GRCm39) annotation onto it with Liftoff.
2. **Nanopore barcode extraction:** builds a splitcode config from the [seqspec](https://github.com/pachterlab/seqspec) in `workflow/configs/`, extracts barcodes, UMI and cDNA from both read orientations, then corrects barcodes against the Parse onlists and merges the oligo-dT and random-hexamer round-1 barcodes from the same well. Each nanopore FASTQ is processed separately; sublibraries are never merged.
3. **Short-read splice junctions:** aligns the Illumina R1 (cDNA) reads with STAR to collect splice junctions.
4. **Transcriptome:** aligns the corrected nanopore cDNA to the genome and builds isoforms with FLAIR 3, using the STAR junctions and the lifted-over annotation.
5. **Quantification:** builds a k=63 kallisto index from the FLAIR isoforms and quantifies each sample per cell with lr-kallisto (`kb count --long`).

## Repository layout

```
8cube_isoform/
├── environment.yaml          # environment for running Snakemake itself
├── workflow/
│   ├── Snakefile
│   ├── envs/                 # one conda environment per tool, used per rule
│   └── configs/
│       ├── nanopore-seqspec.yaml   # read structure of the nanopore libraries
│       ├── config-correct.txt      # splitcode barcode correction
│       ├── config-mergeRT.txt      # splitcode oligo-dT/random-hexamer merge
│       └── onlists/                # Parse barcode onlists
├── resources/                # created by the workflow: genome, annotation, STAR index
├── results/                  # created by the workflow
│   ├── splitcode/<sample>/   # corrected cDNA, UMI and barcode FASTQs
│   ├── star/                 # SJ.out.tab and STAR logs
│   ├── flair/                # aligned BAM and flair.isoforms.{gtf,bed,fa}
│   └── kallisto/             # index and per-sample count matrices
├── logs/                     # one log per rule
└── test/                     # small test dataset and launcher (not tracked)
```

## Running the workflow

### Setup

Create the Snakemake environment once. Each rule's tools are installed automatically from `workflow/envs/`.

```bash
conda env create -f environment.yaml
conda activate 8cube-isoform
```

### Inputs

The input location is set at the top of `workflow/Snakefile` (`RUN_DIR`), or at run time with `--config run_dir=...`. The workflow expects one 8-cube experiment directory laid out as on the lab server:

```
<run_dir>/
├── nanopore/*.fastq.gz                 # one sample per file
└── nova*/<lib>L001_R1_001.fastq.gz     # Illumina R1 (cDNA) reads used for junctions
```

Only lane 1 of the short reads is used, to keep the STAR runtime manageable.

### Run

Run from the repository root; Snakemake finds `workflow/Snakefile` on its own, and loads its default settings (conda environments, cores, memory limit) from `workflow/profiles/default/config.yaml`.

```bash
snakemake
```

Useful options:
- `-n` shows the planned jobs without running anything.
- `--cores N` and `--resources mem_mb=N` override the profile's CPU and memory limits.
- `--config run_dir=/path/to/igvf_0XX` points the workflow at a different experiment.
- `--rerun-incomplete` resumes cleanly after an interruption.

The full run takes a long time (Liftoff and the STAR index alone take hours), so run it inside `tmux` or `screen` when working over SSH.

### Resources

The reference steps are the heaviest: Liftoff uses about 50 GB of RAM on the full GENCODE annotation, and the STAR index about 32 GB. They are built once in `resources/` and reused by later runs.

Each memory-heavy rule declares how much it needs (`resources: mem_mb=...` in the Snakefile), and Snakemake only starts jobs while their total stays under the budget in `workflow/profiles/default/config.yaml` (100 GB by default). **Before running on another machine, set that budget to about 80% of its RAM** (or your allocation on a shared server), and lower `cores` if needed. The workflow needs at least 60 GB for Liftoff.

These limits control scheduling only: they stop Snakemake from starting jobs that wouldn't fit, but don't cap a job that uses more than it declared. For a hard cap on a shared server, also run Snakemake under a memory-limited scope, e.g. `systemd-run --user --scope -p MemoryMax=100G snakemake`.
