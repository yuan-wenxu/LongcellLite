baseCount = function(seq, base = "A") {
  sum(strsplit(seq, "", fixed = TRUE)[[1]] == base)
}

mid_polyA_filter = function(endsites, genome, chr, strand, bin = 20, thresh = 0.4) {
  endsites_uniq = unique(endsites)
  read = tryCatch(
    as.character(get_genome_seq(genome, chr, start = min(endsites_uniq) - bin, end = max(endsites_uniq) + bin, strand = strand)),
    error = function(e) {
      warning(sprintf("Failed to get sequence for '%s': %s", chr, e$message))
      NA_character_
    }
  )
  if (is.na(read)) {
    return(FALSE)
  }

  endsites_offset = endsites_uniq - min(endsites_uniq) + bin
  flag = sapply(endsites_offset, function(i) {
    sub_read = substr(read, start = i - bin + 1, stop = i + bin)
    ratio = baseCount(sub_read, "A") / nchar(sub_read)
    ratio >= thresh
  })
  names(flag) = endsites_uniq
  flag[as.character(endsites)]
}

gene_reads_extraction = function(
  bamFile,
  gene_bed,
  genome,
  toolkit = 5,
  map_qual = 30,
  end_flank = 200,
  splice_site_bin = 2,
  mid_polyA_bin = 20,
  mid_polyA_thresh = 0.4
) {
  chr = unique(gene_bed$chr)[1]
  start = min(gene_bed$start)
  end = max(gene_bed$end)
  exon_bin = as.matrix(gene_bed[, c("start", "end")])
  strand = as.character(unique(gene_bed$strand))

  bam = readBam(bamFile, chr = chr, start = start, end = end, strand = strand, map_qual = map_qual)
  reads = Longcellsrc::extractReads(
    as.character(bam$seq),
    bam$cigar,
    bam$pos,
    exon_bin,
    strand,
    toolkit,
    end_flank,
    splice_site_bin
  )

  reads = as.data.frame(cbind(bam$qname[reads$id], reads %>% dplyr::select(-id)))
  colnames(reads)[1] = "qname"
  reads = reads %>% dplyr::filter(nchar(isoform) > 0)
  if (nrow(reads) == 0) {
    return(reads)
  }

  flags = mid_polyA_filter(reads$isoend, genome, chr, strand, mid_polyA_bin, mid_polyA_thresh)
  reads$polyA[flags] = "0"
  reads
}

reads_extraction = function(
  bam_path,
  gene_bed,
  genome,
  toolkit = 5,
  map_qual = 30,
  end_flank = 200,
  splice_site_bin = 2,
  mid_polyA_bin = 20,
  mid_polyA_thresh = 0.4,
  cores = 1
) {
  genes = unique(gene_bed$gene)
  bamFile = Rsamtools::BamFile(bam_path)
  worker = function(gene_id) {
    sub_bed = gene_bed %>% dplyr::filter(gene == gene_id)
    sub_reads = gene_reads_extraction(
      bamFile = bamFile,
      gene_bed = sub_bed,
      genome = genome,
      toolkit = toolkit,
      map_qual = map_qual,
      end_flank = end_flank,
      splice_site_bin = splice_site_bin,
      mid_polyA_bin = mid_polyA_bin,
      mid_polyA_thresh = mid_polyA_thresh
    )
    if (nrow(sub_reads) == 0) {
      return(NULL)
    }
    sub_reads$gene = gene_id
    sub_reads
  }

  if (cores > 1) {
    reads = parallel::mclapply(genes, worker, mc.cores = cores)
  } else {
    reads = lapply(genes, worker)
  }
  as.data.frame(do.call(rbind, reads))
}

extract_read_isoforms_from_bam = function(
  bam_path,
  gtf_path = NULL,
  gene_bed_path = NULL,
  genome_path,
  genome_name = NULL,
  qname_barcode_pattern = "(?<=CB:)[ACGTN]+",
  qname_umi_pattern = "(?<=UB:)[ACGTN]+",
  toolkit = 5,
  work_dir = "./",
  map_qual = 30,
  end_flank = 200,
  splice_site_bin = 2,
  bedtools = "bedtools",
  samtools = "samtools",
  samtools_threads = 1,
  cores = 1,
  overwrite = FALSE,
  gene_bed = NULL,
  genome = NULL,
  prepared_bam_path = NULL,
  init_work_dir = TRUE
) {
  out_path = file.path(work_dir, "read_isoform", "read_isoform.tsv")
  if (file.exists(out_path) && !overwrite) {
    warning("The read-level isoform output already exists; reusing cached result.")
    return(read.table(out_path, header = TRUE, sep = "\t"))
  }

  if (init_work_dir) {
    init_project(work_dir)
  }
  if (is.null(gene_bed)) {
    annot = annotation(gtf_path = gtf_path, gene_bed_path = gene_bed_path, work_dir = work_dir, overwrite = overwrite)
    gene_bed = annot[[1]]
  }
  if (is.null(genome)) {
    genome = load_genome(genome_name = genome_name, genome_path = genome_path)
  }
  if (is.null(prepared_bam_path)) {
    prepared_bam_path = prepare_input_bam(
      bam_path,
      work_dir = work_dir,
      samtools = samtools,
      samtools_threads = samtools_threads,
      force = overwrite
    )
  }

  gene_range = gene_bed %>% dplyr::group_by(gene) %>%
    dplyr::summarise(chr = unique(chr), start = min(start), end = max(end), strand = unique(strand), .groups = "drop")
  write.table(
    gene_range[, c("chr", "start", "end", "strand", "gene")],
    file.path(work_dir, "annotation", "gene_range.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  noncover = bamGeneCoverage(
    bam = prepared_bam_path,
    gene_range_bed = file.path(work_dir, "annotation", "gene_range.txt"),
    outdir = file.path(work_dir, "annotation"),
    bedtools = bedtools
  )
  if (!is.null(noncover)) {
    gene_bed = gene_bed %>% dplyr::filter(!gene %in% noncover$gene)
  }

  reads = reads_extraction(
    bam_path = prepared_bam_path,
    gene_bed = gene_bed,
    genome = genome,
    toolkit = toolkit,
    map_qual = map_qual,
    end_flank = end_flank,
    splice_site_bin = splice_site_bin,
    cores = cores
  )
  bc = extract_bc_from_qname_vector(
    reads$qname,
    barcode_pattern = qname_barcode_pattern,
    umi_pattern = qname_umi_pattern
  )
  reads_bc = dplyr::inner_join(bc, reads, by = c("name" = "qname")) %>%
    dplyr::mutate(polyA.x = as.numeric(polyA.x), polyA.y = as.numeric(polyA.y)) %>%
    dplyr::mutate(polyA = polyA.x & polyA.y) %>%
    dplyr::select(-polyA.x, -polyA.y)

  saveResult(reads_bc, out_path)
  reads_bc
}
