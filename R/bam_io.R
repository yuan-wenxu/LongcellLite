bamGeneCoverage = function(bam, gene_range_bed, outdir, bedtools = "bedtools") {
  command1 = paste(
    c(bedtools, "bamtobed -i", bam, "|", bedtools, "merge -i - >", file.path(outdir, "cover.bed")),
    collapse = " "
  )
  command2 = paste(
    c(bedtools, "subtract -a", gene_range_bed, "-b", file.path(outdir, "cover.bed"), "-A >", file.path(outdir, "noncover.bed")),
    collapse = " "
  )
  system(command1)
  system(command2)

  filename = file.path(outdir, "noncover.bed")
  if (file.size(filename) == 0L) {
    return(NULL)
  }

  noncover = read.table(filename)
  colnames(noncover) = c("chr", "start", "end", "strand", "gene")
  noncover
}

extract_qname_annotation = function(qnames, pattern, field_name) {
  out = stringr::str_extract(qnames, pattern)
  if (anyNA(out)) {
    missing_n = sum(is.na(out))
    stop(sprintf("Failed to extract %s from %d read names using pattern: %s", field_name, missing_n, pattern))
  }
  out
}

extract_bc_from_qname_vector = function(
  qnames,
  barcode_pattern = "(?<=CB:)[ACGTN]+",
  umi_pattern = "(?<=UB:)[ACGTN]+"
) {
  data.frame(
    name = qnames,
    barcode = extract_qname_annotation(qnames, barcode_pattern, "barcode"),
    umi = extract_qname_annotation(qnames, umi_pattern, "UMI"),
    polyA = 1,
    stringsAsFactors = FALSE
  )
}

prepare_input_bam = function(input_bam_path, work_dir, samtools = "samtools", force = FALSE) {
  if (!file.exists(input_bam_path)) {
    stop(sprintf("Input BAM does not exist: %s", input_bam_path))
  }

  input_bam_path = normalizePath(input_bam_path)
  input_bai_candidates = c(paste0(input_bam_path, ".bai"), sub("\\.bam$", ".bai", input_bam_path))
  has_index = any(file.exists(input_bai_candidates))

  if (!has_index || force) {
    system(paste(c(samtools, "index", input_bam_path), collapse = " "))
  }

  input_bam_path
}

readBam = function(bamFile, chr, start, end, strand, map_qual = 30) {
  if (!strand %in% c("+", "-")) {
    stop("The strand of mapping should be either '+' or '-'!")
  }

  gr = GenomicRanges::GRanges(seqnames = chr, ranges = IRanges::IRanges(start = start, end = end))
  GenomeInfoDb::seqlevels(gr) = GenomeInfoDb::seqlevels(bamFile)

  if (strand == "+") {
    param = Rsamtools::ScanBamParam(
      mapqFilter = map_qual,
      flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE, isMinusStrand = FALSE),
      what = c("qname", "pos", "cigar", "seq"),
      which = gr
    )
  } else {
    param = Rsamtools::ScanBamParam(
      mapqFilter = map_qual,
      flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE, isMinusStrand = TRUE),
      what = c("qname", "pos", "cigar", "seq"),
      which = gr
    )
  }

  aln = Rsamtools::scanBam(bamFile, param = param)
  aln[[1]]
}
