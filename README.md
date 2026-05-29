# LongcellLite

LongcellLite is a lightweight single-cell long-read RNA-seq quantification pipeline built on Longcellsrc for BAM-based isoform and UMI counting.

## Overview

LongcellLite starts from a mapped BAM file with pre-assigned cell barcodes and UMIs stored in the read name. It is designed as a simpler and more maintainable alternative to the heavier `LongcellPre` workflow when the upstream preprocessing and genome alignment have already been completed.

Current scope:

- build gene-level annotation from `GTF` or use an existing gene BED
- read and index a mapped BAM file in place
- reconstruct read-level isoform structures
- extract `CB` and `UB` from read names
- perform UMI clustering and isoform correction
- output a structure-based isoform count table
- optionally map structure-based isoforms to annotated transcripts
- annotate unmatched isoforms as gene-level novel isoforms
- export 10X-style gene and isoform matrices
- write detected-gene and novel-augmented GTF files

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
5. Use the input BAM in place and create an index if needed
6. Filter genes without BAM coverage
7. Extract read-level isoform structures from the BAM
8. Extract barcode and UMI from read names
9. Write `read_isoform/read_isoform.tsv`
10. Perform UMI clustering and isoform correction
11. Write `out/iso_count.tsv`
12. Optionally map structure-based isoforms to annotated transcripts
13. For unmatched isoforms within a known gene, assign novel transcript labels
14. Classify novel transcripts as `nic` or `nnic`
15. Write 10X-style `out/gene/` and `out/isoform/` matrices
16. Write `out/gtf/detected_genes.gtf` and `out/gtf/augmented_with_novel.gtf`

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

SLURM submission is also supported:

```bash
sbatch exec/run_longcelllite.sbatch --config exec/config.sh
```

## Output

Main outputs:

- `annotation/gene_bed.rds`
- `annotation/exon_gtf.rds`
- `read_isoform/read_isoform.tsv`
- `out/iso_count.tsv`

Optional outputs when `--to_isoform` is enabled:

- `out/gene/`
- `out/isoform/`
- `out/gtf/detected_genes.gtf`
- `out/gtf/augmented_with_novel.gtf`

When transcript annotation is enabled, `out/isoform/` may contain:

- known transcript IDs from the input annotation
- novel isoforms labeled as `GeneName.novelX.nic`
- novel isoforms labeled as `GeneName.novelX.nnic`

GTF outputs:

- `out/gtf/detected_genes.gtf`: a filtered version of the input GTF that keeps only genes detected in `out/iso_count.tsv`, plus novel transcript and exon records for those detected genes
- `out/gtf/augmented_with_novel.gtf`: the original input GTF plus transcript and exon records for LongcellLite novel isoforms

Current novel transcript rules:

- `nic`: all exon boundaries can be explained by known exons, but the exon combination is novel
- `nnic`: the isoform contains at least one exon boundary that cannot be explained by the known exon catalog

## Notes

- LongcellLite does not perform adapter detection or barcode extraction from raw FASTQ.
- LongcellLite does not run `minimap2`; it assumes the BAM already exists.
- LongcellLite uses the input BAM in place and does not create a duplicated working-copy BAM.
- The current implementation focuses on the BAM-to-quantification part of the workflow.

## References

- Longcellsrc fork used by this project: [`yuan-wenxu/Longcellsrc`](https://github.com/yuan-wenxu/Longcellsrc)
- Original LongcellPre codebase that inspired this refactor: [`yuntianf/LongcellPre`](https://github.com/yuntianf/LongcellPre)
