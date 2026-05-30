gtf2bed = function(gtf_path, out_path, if_store_binary = TRUE) {
  txdb = txdbmaker::makeTxDbFromGFF(gtf_path, format = "gtf")
  exons_temp = GenomicFeatures::exonicParts(txdb)

  temp = as.data.frame(exons_temp)
  temp = temp[, c(1:5, 8)]
  colnames(temp)[colnames(temp) == "gene_id"] = "gene"
  temp = temp %>% tidyr::unnest(gene)
  temp = as.data.frame(temp)
  temp[, "gene"] = sapply(strsplit(temp[, "gene"], split = ".", fixed = TRUE), function(x) x[1])
  colnames(temp)[1] = "chr"

  temp = temp %>% group_by(gene) %>% mutate(id = exon_id(unique(strand), n()))
  temp = as.data.frame(temp)
  write.table(
    temp,
    file = file.path(out_path, "gene_bed.txt"),
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE
  )
  if (if_store_binary) {
    saveRDS(temp, file = file.path(out_path, "gene_bed.rds"))
  }

  gtf_data = as.data.frame(GenomicFeatures::exons(
    txdb,
    columns = c("gene_id", "tx_name", "exon_name"),
    filter = NULL,
    use.names = FALSE
  ))
  gtf_data = gtf_data %>% tidyr::unnest(gene_id) %>% tidyr::unnest(tx_name)
  colnames(gtf_data)[colnames(gtf_data) == "gene_id"] = "gene"
  gtf_data$gene = sapply(strsplit(gtf_data$gene, split = ".", fixed = TRUE), function(x) x[1])
  gtf_data = gtf_data[, c(2:3, 6:8)]
  colnames(gtf_data)[4:5] = c("transname", "exon_id")

  write.table(
    gtf_data,
    file = file.path(out_path, "exon_gtf.txt"),
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE
  )
  if (if_store_binary) {
    saveRDS(gtf_data, file = file.path(out_path, "exon_gtf.rds"))
  }

  list(temp, gtf_data)
}

exon_id = function(strand, count) {
  if (strand == "+") {
    as.character(seq_len(count))
  } else if (strand == "-") {
    as.character(count:1)
  } else {
    stop("The gene strand should be marked as + or -")
  }
}

createAnnotation = function(
  gtf_path = NULL,
  gene_bed_path = NULL,
  work_dir = "./",
  bed_chr_col = "chr",
  bed_start_col = "start",
  bed_end_col = "end",
  bed_strand_col = "strand",
  bed_gene_col = "gene"
) {
  if (!is.null(gtf_path)) {
    message("Start to build exon annotation based on the gtf file.")
    cache = gtf2bed(gtf_path, file.path(work_dir, "annotation"), if_store_binary = TRUE)
    gene_bed = cache[[1]]
    gtf = cache[[2]]
  } else if (!is.null(gene_bed_path)) {
    gene_bed = read.table(gene_bed_path, header = TRUE)
    gene_bed = gene_bed[, c(bed_chr_col, bed_start_col, bed_end_col, bed_strand_col, bed_gene_col)]
    colnames(gene_bed) = c("chr", "start", "end", "strand", "gene")
    write.table(
      gene_bed,
      file = file.path(work_dir, "annotation", "gene_bed.txt"),
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE
    )
    saveRDS(gene_bed, file = file.path(work_dir, "annotation", "gene_bed.rds"))
    gtf = NULL
  } else {
    stop("Either gtf_path or gene_bed_path must be provided.")
  }

  list(gene_bed, gtf)
}

annotation = function(
  gtf_path = NULL,
  gene_bed_path = NULL,
  work_dir = "./",
  overwrite = FALSE,
  bed_chr_col = "chr",
  bed_start_col = "start",
  bed_end_col = "end",
  bed_strand_col = "strand",
  bed_gene_col = "gene"
) {
  cached_gene_bed = file.path(work_dir, "annotation", "gene_bed.rds")
  cached_gtf = file.path(work_dir, "annotation", "exon_gtf.rds")

  if (file.exists(cached_gene_bed) && !overwrite) {
    cat("Annotation result already exists, skipping this step.\n")
    gene_bed = readRDS(cached_gene_bed)
    gtf = if (file.exists(cached_gtf)) readRDS(cached_gtf) else NULL
    return(list(gene_bed, gtf))
  }

  createAnnotation(
    gtf_path = gtf_path,
    gene_bed_path = gene_bed_path,
    work_dir = work_dir,
    bed_chr_col = bed_chr_col,
    bed_start_col = bed_start_col,
    bed_end_col = bed_end_col,
    bed_strand_col = bed_strand_col,
    bed_gene_col = bed_gene_col
  )
}
