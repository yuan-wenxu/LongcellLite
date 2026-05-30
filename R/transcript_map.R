intron_only = function(reads, gtf, gtf_start_col = "start", gtf_end_col = "end", sep = ",", split = "|") {
  exon_bin = as.data.frame(gtf[, c(gtf_start_col, gtf_end_col)])
  colnames(exon_bin) = c("start", "end")
  exon_bin$chrom = "chr1"

  intron_flag = sapply(reads, function(x) {
    bins = read2bins(x, sep, split)
    bins$chrom = "chr1"
    diff = valr::bed_subtract(bins, exon_bin)
    binsum(diff) == binsum(bins)
  })
  intron_flag
}

binsum = function(bin, start_col = "start", end_col = "end") {
  bin = bin[, c(start_col, end_col)]
  colnames(bin) = c("start", "end")
  bin$chrom = "chr1"
  bin = as.data.frame(valr::bed_merge(bin))
  sum(as.numeric(bin[, "end"]) - as.numeric(bin[, "start"])) + nrow(bin)
}

iso_end_diff <- function(read1, read2, split = "|", sep = ",") {
  read1_bins <- read2bins(read1, split = split, sep = sep)
  read2_bins <- read2bins(read2, split = split, sep = sep)

  start <- max(read1_bins$start[1], read2_bins$start[1])
  end <- min(read1_bins$end[nrow(read1_bins)], read2_bins$end[nrow(read2_bins)])
  left_bound <- min(read1_bins$start[1], read2_bins$start[1])
  right_bound <- max(read1_bins$end[nrow(read1_bins)], read2_bins$end[nrow(read2_bins)])

  if (start == left_bound) {
    left_diff <- 0
  } else {
    left_range <- IRanges::IRanges(start = left_bound, end = start)
    bins1_left <- sum(IRanges::width(IRanges::intersect(IRanges::IRanges(start = read1_bins$start, end = read1_bins$end), left_range)))
    bins2_left <- sum(IRanges::width(IRanges::intersect(IRanges::IRanges(start = read2_bins$start, end = read2_bins$end), left_range)))
    left_diff <- bins1_left - bins2_left
  }

  if (end == right_bound) {
    right_diff <- 0
  } else {
    right_range <- IRanges::IRanges(start = end, end = right_bound)
    bins1_right <- sum(IRanges::width(IRanges::intersect(IRanges::IRanges(start = read1_bins$start, end = read1_bins$end), right_range)))
    bins2_right <- sum(IRanges::width(IRanges::intersect(IRanges::IRanges(start = read2_bins$start, end = read2_bins$end), right_range)))
    right_diff <- bins1_right - bins2_right
  }

  c(left_diff, right_diff)
}

iso_end_diff_v <- Vectorize(iso_end_diff, c("read1", "read2"))

isomatch_penalty = function(left, right) {
  left_penalty = abs(left)
  right_penalty = abs(right)
  if (left > 0) {
    left_penalty = left_penalty^2
  }
  if (right > 0) {
    right_penalty = right_penalty^2
  }
  left_penalty + right_penalty
}

isomatch_penalty_v = Vectorize(isomatch_penalty, c("left", "right"))

penalty2weight = function(x) {
  if (sum(x == 0) > 0) {
    x = ifelse(x > 0, 0, 1)
  } else {
    x = 1 / x
  }
  x / sum(x)
}

nearest_known_site = function(site, known_sites, thresh) {
  if (length(known_sites) == 0) {
    return(list(value = site, known = FALSE))
  }
  diff = abs(known_sites - site)
  idx = which.min(diff)
  if (diff[idx] <= thresh) {
    return(list(value = known_sites[idx], known = TRUE))
  }
  list(value = site, known = FALSE)
}

annotate_novel_isoforms = function(
  transcripts,
  gene,
  gtf,
  thresh = 3,
  gtf_start_col = "start",
  gtf_end_col = "end",
  sep = ",",
  split = "|"
) {
  transcripts_uniq = sort(unique(as.character(transcripts)))
  if (length(transcripts_uniq) == 0) {
    return(NULL)
  }

  known_starts = if (nrow(gtf) > 0) {
    sort(unique(as.numeric(gtf[[gtf_start_col]])))
  } else {
    numeric()
  }
  known_ends = if (nrow(gtf) > 0) {
    sort(unique(as.numeric(gtf[[gtf_end_col]])))
  } else {
    numeric()
  }
  known_exons = if (nrow(gtf) > 0) {
    unique(paste(gtf[[gtf_start_col]], gtf[[gtf_end_col]], sep = sep))
  } else {
    character()
  }

  annotation = lapply(transcripts_uniq, function(tx) {
    bins = read2bins(tx, sep = sep, split = split)
    canon_bins = lapply(seq_len(nrow(bins)), function(i) {
      start_match = nearest_known_site(as.numeric(bins$start[i]), known_starts, thresh)
      end_match = nearest_known_site(as.numeric(bins$end[i]), known_ends, thresh)
      exon_key = paste(start_match$value, end_match$value, sep = sep)
      exon_known = start_match$known && end_match$known && exon_key %in% known_exons
      data.frame(
        start = start_match$value,
        end = end_match$value,
        exon_known = exon_known,
        stringsAsFactors = FALSE
      )
    })
    canon_bins = do.call(rbind, canon_bins)
    class = if (all(canon_bins$exon_known)) "nic" else "nnic"
    signature = paste(paste(canon_bins$start, canon_bins$end, sep = sep), collapse = split)
    data.frame(
      isoform = tx,
      novel_class = class,
      novel_signature = signature,
      stringsAsFactors = FALSE
    )
  })
  annotation = do.call(rbind, annotation)

  novel_groups = unique(annotation[, c("novel_signature", "novel_class")])
  novel_groups = novel_groups[
    order(novel_groups$novel_signature, novel_groups$novel_class),
    ,
    drop = FALSE
  ]
  novel_groups$transname = sprintf(
    "%s.novel%d.%s",
    gene,
    seq_len(nrow(novel_groups)),
    novel_groups$novel_class
  )

  annotation = merge(
    annotation,
    novel_groups,
    by = c("novel_signature", "novel_class"),
    sort = FALSE
  )
  annotation = annotation[
    match(transcripts_uniq, annotation$isoform),
    c("isoform", "transname", "novel_class", "novel_signature")
  ]
  rownames(annotation) = NULL
  annotation
}

novel_index_to_models = function(novel_index, gene, gene_info, split = "|", sep = ",") {
  if (is.null(novel_index) || nrow(novel_index) == 0) {
    return(NULL)
  }
  if (is.null(gene_info) || nrow(gene_info) == 0) {
    return(NULL)
  }

  novel_models = unique(
    novel_index[, c("transname", "novel_class", "novel_signature")]
  )
  novel_models$gene = gene
  novel_models$chr = gene_info$chr[1]
  novel_models$strand = as.character(gene_info$strand[1])
  novel_models$start = vapply(
    strsplit(novel_models$novel_signature, split = split, fixed = TRUE),
    function(x) {
      as.numeric(strsplit(x[1], split = sep, fixed = TRUE)[[1]][1])
    },
    FUN.VALUE = numeric(1)
  )
  novel_models$end = vapply(
    strsplit(novel_models$novel_signature, split = split, fixed = TRUE),
    function(x) {
      last_bin = x[length(x)]
      as.numeric(strsplit(last_bin, split = sep, fixed = TRUE)[[1]][2])
    },
    FUN.VALUE = numeric(1)
  )
  novel_models
}

read_gtf_lines = function(gtf_path) {
  con = if (grepl("\\.gz$", gtf_path, ignore.case = TRUE)) {
    gzfile(gtf_path, open = "rt")
  } else {
    file(gtf_path, open = "rt")
  }
  on.exit(close(con), add = TRUE)
  readLines(con, warn = FALSE)
}

extract_gtf_attribute = function(attr, key) {
  pattern = paste0(key, ' "([^"]+)"')
  match = regexec(pattern, attr)
  values = regmatches(attr, match)
  vapply(values, function(x) {
    if (length(x) >= 2) {
      x[2]
    } else {
      NA_character_
    }
  }, FUN.VALUE = character(1))
}

normalize_gene_id = function(x) {
  sub("\\..*$", "", x)
}

detected_transcript_gtf_lines = function(gtf_path, transcripts, genes = character(), cores = 1) {
  lines = read_gtf_lines(gtf_path)
  header = lines[grepl("^#", lines)]
  body = lines[!grepl("^#", lines)]
  if (length(body) == 0) {
    return(header)
  }

  transcripts = unique(as.character(transcripts))
  transcripts = transcripts[!is.na(transcripts) & nzchar(transcripts)]
  genes = unique(as.character(genes))
  genes = genes[!is.na(genes) & nzchar(genes)]

  chunks = split_indices(length(body), cores)
  chunk_fun = function(idx) {
    sub_body = body[idx]
    fields = strsplit(sub_body, "\t", fixed = TRUE)
    feature = vapply(
      fields,
      function(x) if (length(x) >= 3) x[3] else "",
      FUN.VALUE = character(1)
    )
    attrs = vapply(
      fields,
      function(x) if (length(x) >= 9) x[9] else "",
      FUN.VALUE = character(1)
    )
    gene_ids = normalize_gene_id(extract_gtf_attribute(attrs, "gene_id"))
    transcript_ids = extract_gtf_attribute(attrs, "transcript_id")
    (!is.na(transcript_ids) & transcript_ids %in% transcripts) |
      (feature == "gene" & !is.na(gene_ids) & gene_ids %in% genes)
  }
  if (cores > 1 && length(chunks) > 1) {
    keep = unlist(parallel::mclapply(chunks, chunk_fun, mc.cores = min(cores, length(chunks))), use.names = FALSE)
  } else {
    keep = unlist(lapply(chunks, chunk_fun), use.names = FALSE)
  }

  c(header, body[keep])
}

gtf_gene_metadata = function(gtf_path) {
  lines = read_gtf_lines(gtf_path)
  body = lines[!grepl("^#", lines)]
  if (length(body) == 0) {
    return(data.frame(
      gene = character(),
      gene_name = character(),
      gene_type = character(),
      stringsAsFactors = FALSE
    ))
  }

  fields = strsplit(body, "\t", fixed = TRUE)
  attrs = vapply(
    fields,
    function(x) if (length(x) >= 9) x[9] else "",
    FUN.VALUE = character(1)
  )
  gene = normalize_gene_id(extract_gtf_attribute(attrs, "gene_id"))
  gene_name = extract_gtf_attribute(attrs, "gene_name")
  gene_type = extract_gtf_attribute(attrs, "gene_type")
  gene_type[is.na(gene_type)] = extract_gtf_attribute(
    attrs,
    "gene_biotype"
  )[is.na(gene_type)]

  meta = data.frame(
    gene = gene,
    gene_name = gene_name,
    gene_type = gene_type,
    stringsAsFactors = FALSE
  )
  meta = meta[!is.na(meta$gene) & nzchar(meta$gene), , drop = FALSE]
  meta = meta[!duplicated(meta$gene), , drop = FALSE]
  meta
}

split_gtf_by_gene = function(
  gtf,
  gene_col = "gene",
  iso_col = "transname",
  start_col = "start",
  end_col = "end"
) {
  if (is.null(gtf) || nrow(gtf) == 0) {
    return(list())
  }

  gtf = as.data.frame(gtf) %>%
    dplyr::arrange(.data[[gene_col]], .data[[iso_col]], .data[[start_col]], .data[[end_col]])
  split(gtf, as.character(gtf[[gene_col]]))
}

build_gene_info_index = function(
  gene_bed,
  gene_col = "gene",
  chr_col = "chr",
  strand_col = "strand"
) {
  if (is.null(gene_bed) || nrow(gene_bed) == 0) {
    return(list())
  }

  info = gene_bed %>%
    dplyr::group_by(.data[[gene_col]]) %>%
    dplyr::summarise(
      chr = dplyr::first(.data[[chr_col]]),
      strand = as.character(dplyr::first(.data[[strand_col]])),
      .groups = "drop"
  )
  split(info, as.character(info[[gene_col]]))
}

split_indices = function(n, chunks) {
  if (n <= 0) {
    return(list(integer()))
  }

  chunks = max(1L, min(as.integer(chunks), n))
  if (chunks == 1L) {
    return(list(seq_len(n)))
  }
  split(seq_len(n), cut(seq_len(n), breaks = chunks, labels = FALSE))
}

novel_models_to_gtf_lines = function(
  novel_models,
  gene_meta,
  source = "LongcellLite",
  split = "|",
  sep = ",",
  cores = 1
) {
  if (is.null(novel_models) || nrow(novel_models) == 0) {
    return(character())
  }

  gene_meta = as.data.frame(gene_meta)
  if (nrow(gene_meta) == 0) {
    gene_meta = data.frame(
      gene = character(),
      gene_name = character(),
      gene_type = character(),
      stringsAsFactors = FALSE
    )
  }

  chunks = split_indices(nrow(novel_models), cores)
  chunk_fun = function(idx) {
    unlist(lapply(idx, function(i) {
      row = novel_models[i, , drop = FALSE]
      meta = gene_meta[match(row$gene, gene_meta$gene), , drop = FALSE]
      gene_name = if (nrow(meta) == 1 &&
          !is.na(meta$gene_name) &&
          nzchar(meta$gene_name)) {
        meta$gene_name
      } else {
        row$gene
      }
      gene_type = if (nrow(meta) == 1 &&
          !is.na(meta$gene_type) &&
          nzchar(meta$gene_type)) {
        meta$gene_type
      } else {
        "novel_transcript"
      }

      transcript_attr = paste0(
        'gene_id "', row$gene, '"; ',
        'transcript_id "', row$transname, '"; ',
        'gene_name "', gene_name, '"; ',
        'transcript_name "', row$transname, '"; ',
        'gene_type "', gene_type, '"; ',
        'transcript_type "', row$novel_class, '"; ',
        'novel_class "', row$novel_class, '";'
      )

      transcript_line = paste(
        row$chr, source, "transcript", row$start, row$end, ".", row$strand, ".",
        transcript_attr,
        sep = "\t"
      )

      bins = read2bins(as.character(row$novel_signature), sep = sep, split = split)
      exon_lines = vapply(seq_len(nrow(bins)), function(j) {
        exon_attr = paste0(
          'gene_id "', row$gene, '"; ',
          'transcript_id "', row$transname, '"; ',
          'gene_name "', gene_name, '"; ',
          'transcript_name "', row$transname, '"; ',
          'gene_type "', gene_type, '"; ',
          'transcript_type "', row$novel_class, '"; ',
          'novel_class "', row$novel_class, '"; ',
          'exon_number ', j, '; ',
          'exon_id "', row$transname, '.exon', j, '";'
        )
        paste(
          row$chr, source, "exon", bins$start[j], bins$end[j], ".", row$strand, ".",
          exon_attr,
          sep = "\t"
        )
      }, FUN.VALUE = character(1))

      c(transcript_line, exon_lines)
    }), use.names = FALSE)
  }
  if (cores > 1 && length(chunks) > 1) {
    lines = unlist(parallel::mclapply(chunks, chunk_fun, mc.cores = min(cores, length(chunks))), use.names = FALSE)
  } else {
    lines = unlist(lapply(chunks, chunk_fun), use.names = FALSE)
  }

  lines
}

write_augmented_gtf = function(gtf_path, novel_models, out_path, split = "|", sep = ",", cores = 1) {
  original_lines = read_gtf_lines(gtf_path)
  novel_lines = novel_models_to_gtf_lines(
    novel_models,
    gtf_gene_metadata(gtf_path),
    split = split,
    sep = sep,
    cores = cores
  )
  writeLines(c(original_lines, novel_lines), out_path)
  invisible(NULL)
}

write_detected_transcripts_gtf = function(gtf_path, transcripts, genes, out_path, cores = 1) {
  writeLines(
    detected_transcript_gtf_lines(
      gtf_path,
      transcripts = transcripts,
      genes = genes,
      cores = cores
    ),
    out_path
  )
  invisible(NULL)
}

iso_corres = function(transcripts, gene, gtf, thresh = 3, overlap_thresh = 0,
                      end_bias = 200, gtf_gene_col = "gene", gtf_iso_col = "transname",
                      gtf_start_col = "start", gtf_end_col = "end",
                      sep = ",", split = "|") {
  sub_gtf = gtf %>%
    dplyr::filter(.data[[gtf_gene_col]] == gene) %>%
    dplyr::arrange(.data[[gtf_iso_col]], .data[[gtf_start_col]], .data[[gtf_end_col]])

  sub_gtf_iso = sub_gtf %>%
    dplyr::group_by(.data[[gtf_iso_col]]) %>%
    dplyr::reframe(iso = paste(paste(.data[[gtf_start_col]], .data[[gtf_end_col]], sep = sep), collapse = split))
  sub_gtf_iso = as.data.frame(sub_gtf_iso)
  transcripts_uniq = unique(transcripts)

  transcripts_iso_corres = Longcellsrc::isoset_mid_diff(
    transcripts_uniq, sub_gtf_iso$iso, thresh, overlap_thresh, end_bias, split, sep
  )
  if (transcripts_iso_corres[1, 1] == -1) {
    return(NULL)
  }

  transcripts_iso_corres = as.data.frame(transcripts_iso_corres)
  colnames(transcripts_iso_corres) = c("isoform", "transname", "dis", "overlap")
  suppressWarnings({
    transcripts_iso_corres = transcripts_iso_corres %>%
      dplyr::group_by(isoform) %>%
      dplyr::filter(dis == min(dis)) %>%
      dplyr::select(isoform, transname, overlap) %>%
      dplyr::arrange(isoform, dplyr::desc(overlap))
  })
  transcripts_iso_corres = as.data.frame(transcripts_iso_corres)
  end_diff = as.data.frame(t(iso_end_diff_v(
    transcripts_uniq[transcripts_iso_corres$isoform + 1],
    sub_gtf_iso[transcripts_iso_corres$transname + 1, "iso"]
  )))
  colnames(end_diff) = c("left", "right")
  rownames(end_diff) = NULL

  transcripts_iso_corres = cbind(transcripts_iso_corres, end_diff)
  transcripts_iso_corres = transcripts_iso_corres %>%
    dplyr::mutate(
      isoform = transcripts_uniq[isoform + 1],
      transname = sub_gtf_iso[transname + 1, gtf_iso_col]
    )

  transcripts_iso_corres = transcripts_iso_corres %>%
    dplyr::mutate(penalty = isomatch_penalty_v(left, right)) %>%
    dplyr::group_by(isoform) %>%
    dplyr::mutate(levels = length(unique(transname))) %>%
    dplyr::arrange(levels, penalty)

  transcripts_iso_corres_filter = transcripts_iso_corres %>% dplyr::filter(levels < 3)
  if (nrow(transcripts_iso_corres_filter) > 0) {
    isos = unique(transcripts_iso_corres_filter$transname)
    transcripts_iso_corres = transcripts_iso_corres %>%
      dplyr::filter(transname %in% isos) %>%
      dplyr::group_by(isoform) %>%
      dplyr::mutate(levels = length(unique(transname)))
  }
  transcripts_iso_corres %>% dplyr::group_by(isoform) %>% dplyr::mutate(weight = penalty2weight(penalty))
}

cells_genes_isos_count = function(data, gtf, thresh = 3, overlap_thresh = 0,
                                  filter_only_intron = TRUE, cell_col = "cell", gene_col = "gene",
                                  transcript_col = "isoform", count_col = "count",
                                  gtf_gene_col = "gene", gtf_iso_col = "transname",
                                  gtf_start_col = "start", gtf_end_col = "end",
                                  split = "|", sep = ",",
                                  data_by_gene = NULL, gtf_by_gene = NULL,
                                  gene_info_by_gene = NULL, return_novel_models = FALSE) {
  data = as.data.frame(data)
  gene_uniq = unique(data[, gene_col])
  if (is.null(data_by_gene)) {
    data_by_gene = split(data, as.character(data[[gene_col]]))
  }
  if (is.null(gtf_by_gene)) {
    gtf_by_gene = split_gtf_by_gene(
      gtf,
      gene_col = gtf_gene_col,
      iso_col = gtf_iso_col,
      start_col = gtf_start_col,
      end_col = gtf_end_col
    )
  }

  out = lapply(gene_uniq, function(i) {
    sub_data = data_by_gene[[i]]
    sub_gtf = gtf_by_gene[[i]]
    if (is.null(sub_gtf)) {
      sub_gtf = gtf[0, , drop = FALSE]
    }
    gene_info = NULL
    if (!is.null(gene_info_by_gene)) {
      gene_info = gene_info_by_gene[[i]]
    }

    if (filter_only_intron && nrow(sub_gtf) > 0) {
      transcripts_uniq = unique(sub_data[, transcript_col])
      intron_flag = intron_only(transcripts_uniq, sub_gtf, gtf_start_col, gtf_end_col, sep, split)
      transcripts_uniq = transcripts_uniq[!intron_flag]
      if (length(transcripts_uniq) == 0) {
        return(NULL)
      }
      sub_data = sub_data %>% dplyr::filter(.data[[transcript_col]] %in% transcripts_uniq)
    }

    iso_index = NULL
    if (nrow(sub_gtf) > 0) {
      iso_index = iso_corres(
        sub_data[, transcript_col], gene = i, gtf = sub_gtf, thresh = thresh,
        overlap_thresh = overlap_thresh, gtf_gene_col = gtf_gene_col,
        gtf_iso_col = gtf_iso_col, gtf_start_col = gtf_start_col,
        gtf_end_col = gtf_end_col, sep = sep, split = split
      )
    }

    matched_isoforms = if (is.null(iso_index)) character() else unique(as.character(iso_index$isoform))
    known_sub_data = sub_data %>% dplyr::filter(.data[[transcript_col]] %in% matched_isoforms)
    novel_sub_data = sub_data %>% dplyr::filter(!(.data[[transcript_col]] %in% matched_isoforms))

    known_out = NULL
    novel_models = NULL
    if (!is.null(iso_index) && nrow(known_sub_data) > 0) {
      iso_index_weight = tidyr::pivot_wider(
        iso_index[, c("isoform", "weight", "transname")],
        names_from = "transname",
        values_from = "weight"
      )
      iso_index_weight[is.na(iso_index_weight)] = 0

      known_out = cbind(
        known_sub_data$cell,
        iso_index_weight[match(known_sub_data$isoform, iso_index_weight$isoform), 2:ncol(iso_index_weight)] * known_sub_data$count
      )
      colnames(known_out)[1] = "cell"
      known_out = tidyr::pivot_longer(
        known_out,
        cols = setdiff(colnames(known_out), "cell"),
        names_to = "isoform",
        values_to = "count"
      )
      known_out = known_out %>% dplyr::filter(count > 0)
    }

    novel_out = NULL
    if (nrow(novel_sub_data) > 0) {
      novel_index = annotate_novel_isoforms(
        novel_sub_data[, transcript_col],
        gene = i,
        gtf = sub_gtf,
        thresh = thresh,
        gtf_start_col = gtf_start_col,
        gtf_end_col = gtf_end_col,
        sep = sep,
        split = split
      )
      novel_out = dplyr::left_join(
        novel_sub_data[, c(cell_col, transcript_col, count_col)],
        novel_index,
        by = setNames("isoform", transcript_col)
      ) %>%
        dplyr::group_by(.data[[cell_col]], transname) %>%
        dplyr::summarise(count = sum(.data[[count_col]]), .groups = "drop")
      colnames(novel_out) = c("cell", "isoform", "count")
      novel_models = novel_index_to_models(novel_index, i, gene_info, split = split, sep = sep)
    }

    sub_out = dplyr::bind_rows(known_out, novel_out)
    if (is.null(sub_out) || nrow(sub_out) == 0) {
      return(list(count = NULL, novel_models = novel_models))
    }
    sub_out$gene = i
    list(count = sub_out, novel_models = novel_models)
  })
  if (return_novel_models) {
    count_out = Filter(Negate(is.null), lapply(out, `[[`, "count"))
    novel_out = Filter(Negate(is.null), lapply(out, `[[`, "novel_models"))
    return(list(
      count = if (length(count_out) > 0) as.data.frame(do.call(dplyr::bind_rows, count_out)) else NULL,
      novel_models = if (length(novel_out) > 0) {
        as.data.frame(do.call(dplyr::bind_rows, novel_out)) %>% dplyr::distinct()
      } else {
        NULL
      }
    ))
  }
  do.call(rbind, lapply(out, `[[`, "count"))
}

UMI_count_to_isoform = function(umi_count, dir, gene_bed, gtf = NULL, gene_col = "gene",
                                bed_gene_col = "gene", bed_strand_col = "strand",
                                filter_only_intron = TRUE, mid_offset_thresh = 3,
                                overlap_thresh = 0, gtf_gene_col = "gene",
                                gtf_start_col = "start", gtf_end_col = "end",
                                gtf_iso_col = "transname", split = "|", sep = ",",
                                cores = 1, gtf_source_path = NULL) {
  cat("Start to do isoform alignment:\n")
  if (is.null(gtf)) {
    warning("The gtf annotation is not provided for the isoform imputation, will skip this step!")
    return(NULL)
  }

  cores = coreDetect(cores)
  novel_models = NULL
  count_mat_rds = file.path(dir, "isoform_count_mat.rds")
  gtf_dir = file.path(dir, "gtf")
  novel_models_rds = file.path(gtf_dir, "novel_transcript_models.rds")
  gtf_by_gene = split_gtf_by_gene(
    gtf,
    gene_col = gtf_gene_col,
    iso_col = gtf_iso_col,
    start_col = gtf_start_col,
    end_col = gtf_end_col
  )
  gene_info_by_gene = build_gene_info_index(
    gene_bed,
    gene_col = bed_gene_col,
    chr_col = "chr",
    strand_col = bed_strand_col
  )
  if (file.exists(count_mat_rds) && file.exists(novel_models_rds)) {
    cat("Isoform matrix cache and novel transcript models already exist, skipping this step.\n")
    count_mat = readRDS(count_mat_rds)
    novel_models = readRDS(novel_models_rds)
  } else {
    data_split = genes_distribute(umi_count, 16 * cores, gene_col)
    count_mat_fun = function(x) {
      if (is.null(x)) {
        return(NULL)
      }
      data_by_gene = split(as.data.frame(x), as.character(x[[gene_col]]))
      sub_result = cells_genes_isos_count(
        x, gtf,
        thresh = mid_offset_thresh,
        overlap_thresh = overlap_thresh,
        filter_only_intron = filter_only_intron,
        gtf_gene_col = gtf_gene_col,
        gtf_iso_col = gtf_iso_col,
        gtf_start_col = gtf_start_col,
        gtf_end_col = gtf_end_col,
        split = split, sep = sep,
        data_by_gene = data_by_gene,
        gtf_by_gene = gtf_by_gene,
        gene_info_by_gene = gene_info_by_gene,
        return_novel_models = TRUE
      )
      if (is.null(sub_result$count) || nrow(sub_result$count) == 0) {
        return(NULL)
      }
      sub_result
    }

    if (cores > 1) {
      count_mat = parallel::mclapply(data_split, count_mat_fun, mc.cores = cores)
    } else {
      count_mat = lapply(data_split, count_mat_fun)
    }
    count_out = Filter(Negate(is.null), lapply(count_mat, `[[`, "count"))
    novel_out = Filter(Negate(is.null), lapply(count_mat, `[[`, "novel_models"))
    count_mat = as.data.frame(do.call(dplyr::bind_rows, count_out))
    novel_models = if (length(novel_out) > 0) {
      as.data.frame(do.call(dplyr::bind_rows, novel_out)) %>% dplyr::distinct()
    } else {
      NULL
    }
    saveIsoMat(count_mat, dir)
    saveRDS(count_mat, count_mat_rds)
  }

  if (!is.null(gtf_source_path) && nzchar(gtf_source_path)) {
    dir.create(gtf_dir, showWarnings = FALSE, recursive = TRUE)

    if (!is.null(novel_models)) {
      saveRDS(novel_models, novel_models_rds)
    }

    write_augmented_gtf(
      gtf_path = gtf_source_path,
      novel_models = novel_models,
      out_path = file.path(gtf_dir, "augmented_with_novel.gtf"),
      split = split,
      sep = sep,
      cores = cores
    )

    expressed_transcripts = unique(as.character(count_mat$isoform))
    expressed_genes = unique(as.character(count_mat$gene))
    write_detected_transcripts_gtf(
      gtf_path = file.path(gtf_dir, "augmented_with_novel.gtf"),
      transcripts = expressed_transcripts,
      genes = expressed_genes,
      out_path = file.path(gtf_dir, "detected_genes.gtf"),
      cores = cores
    )
  }

  attr(count_mat, "novel_models") = novel_models
  attr(count_mat, "structure_output") = novel_models
  invisible(count_mat)
}
