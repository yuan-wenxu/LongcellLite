#!/usr/bin/env Rscript

library(LongcellLite)
library(argparse)

p = ArgumentParser(description = "LongcellLite: BAM-based isoform and UMI quantification")
p$add_argument("--bam", required = TRUE, help = "Input mapped BAM file")
p$add_argument("--gtf", required = FALSE, help = "Input GTF annotation")
p$add_argument("--gene_bed_path", required = FALSE, help = "Optional prebuilt gene BED annotation")
p$add_argument("--genome_path", required = TRUE, help = "Reference genome FASTA")
p$add_argument("--genome_name", required = FALSE, help = "Optional BSgenome-compatible genome name")
p$add_argument("--qname_barcode_pattern", default = "(?<=CB:)[ACGTN]+", help = "Regex for barcode in qname")
p$add_argument("--qname_umi_pattern", default = "(?<=UB:)[ACGTN]+", help = "Regex for UMI in qname")
p$add_argument("--toolkit", default = 5, type = "integer", help = "Library orientation, 5 or 3")
p$add_argument("--work_dir", default = "./", help = "Output directory")
p$add_argument("--map_qual", default = 30, type = "integer", help = "Minimum mapping quality")
p$add_argument("--end_flank", default = 200, type = "integer", help = "Flanking region size for TSS/TES assignment")
p$add_argument("--splice_site_bin", default = 2, type = "integer", help = "Bin size for splice site assignment")
p$add_argument("--samtools", default = "samtools", help = "samtools binary")
p$add_argument("--bedtools", default = "bedtools", help = "bedtools binary")
p$add_argument("--cores", default = 1, type = "integer", help = "Number of cores")
p$add_argument("--overwrite", action = "store_true", help = "Overwrite cached outputs")
p$add_argument("--to_isoform", action = "store_true", help = "Generate gene/isoform matrices from GTF mapping")
p$add_argument("--splice_site_thresh", default = 3, type = "integer", help = "Threshold for splice site assignment")
p$add_argument("--mid_offset_thresh", default = 5, type = "integer", help = "Threshold for mid offset assignment")
p$add_argument("--overlap_thresh", default = 0, type = "double", help = "Threshold for overlap assignment")
p$add_argument("--filter_only_intron", action = "store_true", help = "Filter only intron reads")

args = p$parse_args()

run_longcelllite(
  bam_path = args$bam,
  gtf_path = args$gtf,
  gene_bed_path = args$gene_bed_path,
  genome_path = args$genome_path,
  genome_name = args$genome_name,
  qname_barcode_pattern = args$qname_barcode_pattern,
  qname_umi_pattern = args$qname_umi_pattern,
  toolkit = args$toolkit,
  work_dir = args$work_dir,
  map_qual = args$map_qual,
  end_flank = args$end_flank,
  splice_site_bin = args$splice_site_bin,
  bedtools = args$bedtools,
  samtools = args$samtools,
  cores = args$cores,
  overwrite = isTRUE(args$overwrite),
  to_isoform = isTRUE(args$to_isoform),
  splice_site_thresh = args$splice_site_thresh,
  mid_offset_thresh = args$mid_offset_thresh,
  overlap_thresh = args$overlap_thresh,
  filter_only_intron = isTRUE(args$filter_only_intron)
)
