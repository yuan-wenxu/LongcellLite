isoform_correct_filter <- function(gene_cells_cluster, filter_ratio = 0, strand, split = "|", sep = ",") {
  gene_isoform = splice_site_table(
    gene_cells_cluster$isoform,
    split, sep,
    splice_site_thresh = 0
  )
  if (length(gene_isoform) == 0 || nrow(gene_isoform) == 0) {
    return(NULL)
  }

  gene_isoform = gene_isoform %>% dplyr::select(-id)
  if (ncol(gene_isoform) > 2) {
    filter = as.data.frame(gene_isoform %>% dplyr::select(-c(start, end)))
    filter = colSums(filter, na.rm = TRUE) > 0
    gene_isoform = gene_isoform[, c(TRUE, filter, TRUE)]
  }

  gene_isoform = cells_isoform_correct(
    gene_cells_cluster$cell,
    gene_cells_cluster$cluster,
    gene_isoform,
    gene_cells_cluster$polyA
  )
  if (nrow(gene_isoform) == 0) {
    return(gene_isoform)
  }
  as.data.frame(cells_isoforms_size_filter(gene_isoform, ratio = filter_ratio))
}

gene_umi_count <- function(cell_exon, strand, bar = "barcode",
                           isoform = "isoform", polyA = "polyA",
                           sim_thresh = NULL, split = "|", sep = ",",
                           splice_site_thresh = 3, verbose = FALSE,
                           filter_ratio = 0) {
  colnames(cell_exon)[which(colnames(cell_exon) == bar)] = "cell"
  colnames(cell_exon)[which(colnames(cell_exon) == isoform)] = "isoform"
  colnames(cell_exon)[which(colnames(cell_exon) == polyA)] = "polyA"
  cell_exon$polyA = as.numeric(cell_exon$polyA)

  gene_isoform = splice_site_table(cell_exon$isoform, split, sep, splice_site_thresh)
  if (length(gene_isoform) == 0 || nrow(gene_isoform) == 0) {
    return(NULL)
  }
  cell_exon = cell_exon[gene_isoform$id, ]

  cells = unique(cell_exon[, "cell"])
  if (is.null(sim_thresh)) {
    sim_thresh = nchar(cell_exon$umi)[1] / 2 + 1
  }

  gene_cells_cluster <- lapply(cells, function(i) {
    cell_i = cell_exon[cell_exon[, "cell"] == i, ]
    if (verbose) {
      cat(nrow(cell_i), " reads in cell ", i, "\n")
    }

    cell_i$cluster = 0
    if (nrow(cell_i) > 400000 || length(unique(cell_i$umi)) > 20000) {
      warning("Too many reads in cell ", i, " which exceeds the max limit of memory")
      return(NULL)
    }
    if (nrow(cell_i) != 1) {
      cell_i$cluster = umi_cluster(cell_i$umi, iso = cell_i$isoform, thresh = sim_thresh)
    } else {
      cell_i$cluster = 1
    }
    cell_i
  })

  gene_cells_cluster = as.data.frame(do.call(rbind, gene_cells_cluster))

  isoform_correct_filter(gene_cells_cluster, filter_ratio, strand, split = split, sep = sep)
}

umi_count <- function(cell_exon, gene_strand, bar = "barcode", gene = "gene",
                      isoform = "isoform", polyA = "polyA", sim_thresh = NULL,
                      split = "|", sep = ",", splice_site_thresh = 3, verbose = FALSE,
                      filter_ratio = 0) {
  genes <- unique(as.character(cell_exon[[gene]]))
  genes_umi_count <- lapply(genes, function(i) {
    if (verbose) {
      cat(i, "\n")
    }
    sub_cell_exon = cell_exon[cell_exon[, gene] == i, ]
    strand = unique(gene_strand[gene_strand$gene == i, "strand"])
    if (nrow(sub_cell_exon) < splice_site_thresh) {
      if (verbose) {
        cat("too few reads, will be filtered out\n")
      }
      return(NULL)
    }
    tryCatch({
      sub_umi_count = gene_umi_count(
        sub_cell_exon, strand = strand,
        bar = bar, isoform = isoform, polyA = polyA,
        sim_thresh = sim_thresh, split = split, sep = sep,
        splice_site_thresh = splice_site_thresh, verbose = verbose,
        filter_ratio = filter_ratio
      )
      if (is.null(sub_umi_count) || nrow(sub_umi_count) == 0) {
        return(NULL)
      }
      sub_umi_count$gene = i
      sub_umi_count
    }, error = function(e) {
      message("Error processing gene ", i, ": ", conditionMessage(e))
      NULL
    })
  })

  do.call(rbind, genes_umi_count)
}

umi_count_parallel <- function(data, dir, gene_bed,
                               bar = "barcode", gene = "gene",
                               isoform = "isoform", polyA = "polyA",
                               sim_thresh = NULL, split = "|", sep = ",",
                               splice_site_thresh = 3, verbose = FALSE,
                               bed_gene_col = "gene", bed_strand_col = "strand",
                               cores = 1, force_UMI_dedup = FALSE,
                               filter_ratio = 0) {
  cat("Start to do UMI deduplication:\n")

  out_path = file.path(dir, "iso_count.tsv")
  if (file.exists(out_path) && !force_UMI_dedup) {
    warning("The UMI deduplication output already exists; reusing cached result.")
    return(read.table(out_path, header = TRUE, sep = "\t"))
  }

  cores = coreDetect(cores)
  data_split = genes_distribute(data, 16 * cores, gene)
  gene_strand = unique(gene_bed[, c(bed_gene_col, bed_strand_col)])

  count_fun = function(x) {
    sub_count = umi_count(
      x, gene_strand,
      bar = bar, gene = gene, isoform = isoform, polyA = polyA,
      sim_thresh = sim_thresh, split = split, sep = sep,
      splice_site_thresh = splice_site_thresh, verbose = verbose,
      filter_ratio = filter_ratio
    )
    if (length(sub_count) == 0 || nrow(sub_count) == 0) {
      return(NULL)
    }
    sub_count
  }

  if (cores > 1) {
    count = parallel::mclapply(data_split, count_fun, mc.cores = cores)
  } else {
    count = lapply(data_split, count_fun)
  }

  count = as.data.frame(do.call(rbind, count))
  count = count %>% dplyr::select(cell, gene, isoform, count, polyA)
  saveResult(count, out_path)
  count
}

run_longcelllite <- function(bam_path,
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
                             cores = 1,
                             overwrite = FALSE,
                             to_isoform = TRUE,
                             splice_site_thresh = 3,
                             mid_offset_thresh = 3,
                             overlap_thresh = 0,
                             filter_only_intron = TRUE) {
  init_project(work_dir)
  annot = annotation(gtf_path = gtf_path, gene_bed_path = gene_bed_path, work_dir = work_dir, overwrite = overwrite)
  gene_bed = annot[[1]]
  gtf = annot[[2]]
  genome = load_genome(genome_name = genome_name, genome_path = genome_path)
  prepared_bam_path = prepare_input_bam(bam_path, work_dir = work_dir, samtools = samtools, force = overwrite)

  reads_bc = extract_read_isoforms_from_bam(
    bam_path = bam_path,
    gtf_path = gtf_path,
    gene_bed_path = gene_bed_path,
    genome_path = genome_path,
    genome_name = genome_name,
    qname_barcode_pattern = qname_barcode_pattern,
    qname_umi_pattern = qname_umi_pattern,
    toolkit = toolkit,
    work_dir = work_dir,
    map_qual = map_qual,
    end_flank = end_flank,
    splice_site_bin = splice_site_bin,
    bedtools = bedtools,
    samtools = samtools,
    cores = cores,
    overwrite = overwrite,
    gene_bed = gene_bed,
    genome = genome,
    prepared_bam_path = prepared_bam_path,
    init_work_dir = FALSE
  )

  iso_count = umi_count_parallel(
    data = reads_bc,
    dir = file.path(work_dir, "out"),
    gene_bed = gene_bed,
    cores = cores,
    splice_site_thresh = splice_site_thresh,
    force_UMI_dedup = overwrite
  )

  if (to_isoform) {
    UMI_count_to_isoform(
      umi_count = iso_count,
      dir = file.path(work_dir, "out"),
      gene_bed = gene_bed,
      gtf = gtf,
      gtf_source_path = gtf_path,
      cores = cores,
      filter_only_intron = filter_only_intron,
      mid_offset_thresh = mid_offset_thresh,
      overlap_thresh = overlap_thresh
    )
  }

  invisible(list(read_isoform = reads_bc, iso_count = iso_count))
}
