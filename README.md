# LongcellLite

LongcellLite is a lightweight single-cell long-read RNA-seq quantification pipeline built on Longcellsrc for BAM-based isoform and UMI counting.

## Overview

LongcellLite starts from a mapped BAM file with pre-assigned cell barcodes and UMIs stored in the read name. It is designed as a simpler and more maintainable alternative to the heavier `LongcellPre` workflow when the upstream preprocessing and genome alignment have already been completed.

Current scope:

- build gene-level annotation from `GTF` or use an existing gene BED
- read a mapped BAM file
- reconstruct read-level isoform structures
- extract `CB` and `UB` from read names
- perform UMI clustering and isoform correction
- output a structure-based isoform count table
- optionally map structure-based isoforms to annotated transcripts and export 10X-style matrices

## Input

LongcellLite expects:

- a mapped BAM file
- a reference genome FASTA
- a GTF annotation or a prebuilt gene BED
- read names containing cell barcodes and UMIs, for example:

```text
read123 CB:ACGTACGTACGTACGT UB:TTGCAAAAAAAA
```

Default patterns:

- barcode: `(?<=CB:)[ACGTN]+`
- UMI: `(?<=UB:)[ACGTN]+`

## Workflow

The current first-version workflow is:

1. Parse CLI arguments from `exec/RunLongcellLite.R`
2. Initialize the output directory structure
3. Build or load annotation
4. Load the reference genome
5. Prepare and index the input BAM
6. Filter genes without BAM coverage
7. Extract read-level isoform structures from the BAM
8. Extract barcode and UMI from read names
9. Write `read_isoform/read_isoform.tsv`
10. Perform UMI clustering and isoform correction
11. Write `out/iso_count.tsv`
12. Optionally map structure-based isoforms to annotated transcripts
13. Optionally write `out/gene/` and `out/isoform/` 10X-style matrices

## Installation

### 1. Create the pixi environment

```bash
cd ./LongcellLite
pixi install
```

### 2. Fix `GenomeInfoDbData` in the pixi R environment

Some `bioconductor-genomeinfodbdata` builds install incomplete metadata through the post-link step. LongcellLite provides a repair task for that package:

```bash
pixi run install-genomeinfodbdata
```

### 3. Install `Longcellsrc`

LongcellLite reuses `Longcellsrc` for core low-level computations:

```bash
pixi run install-longcellsrc
```

This currently installs:

- [`yuan-wenxu/Longcellsrc`](https://github.com/yuan-wenxu/Longcellsrc)


### 4. Install LongcellLite itself

```bash
pixi run install-longcelllite-local
```

### 5. Install `Longcellsrc` and LongcellLite itself

```bash
pixi run install-longcelllite-stack
```

## Usage

Example:

```bash
pixi run -e long Rscript exec/RunLongcellLite.R \
  --bam /path/to/mapping.bam \
  --gtf /path/to/annotation.gtf.gz \
  --genome_path /path/to/genome.fa.gz \
  --work_dir /path/to/output \
  --cores 8 \
  --to_isoform
```

## Output

Main outputs:

- `annotation/gene_bed.rds`
- `annotation/exon_gtf.rds`
- `bam/polish.bam`
- `read_isoform/read_isoform.tsv`
- `out/iso_count.tsv`

Optional outputs when `--to_isoform` is enabled:

- `out/gene/`
- `out/isoform/`

## Notes

- LongcellLite does not perform adapter detection or barcode extraction from raw FASTQ.
- LongcellLite does not run `minimap2`; it assumes the BAM already exists.
- The current implementation focuses on the BAM-to-quantification part of the workflow.

## References

- Longcellsrc fork used by this project: [`yuan-wenxu/Longcellsrc`](https://github.com/yuan-wenxu/Longcellsrc)
- Original LongcellPre codebase that inspired this refactor: [`yuntianf/LongcellPre`](https://github.com/yuntianf/LongcellPre)
