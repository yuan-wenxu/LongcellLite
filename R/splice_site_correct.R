splice_site_table <- function(isoforms, split = "|", sep = ",", splice_site_thresh = 10) {
  out <- Longcellsrc::splice_site_table_cpp(isoforms, split, sep, splice_site_thresh)
  out$start <- as.numeric(out$start)
  out$end <- as.numeric(out$end)
  out = as.data.frame(do.call(cbind, out))

  if (ncol(out) > 3) {
    mid_names <- setdiff(colnames(out), c("id", "start", "end"))
    ord <- order(suppressWarnings(as.numeric(mid_names)))
    mid_cols <- mid_names[ord]

    lower <- out$start
    upper <- out$end

    for (nm in mid_cols) {
      v <- suppressWarnings(as.numeric(out[[nm]]))
      bad <- (as.numeric(nm) < lower) | (as.numeric(nm) > upper)
      if (any(bad, na.rm = TRUE)) {
        v[bad] <- NA_real_
      }
      out[[nm]] <- v
    }
    out <- out[c("id", "start", mid_cols, "end")]
  }

  rownames(out) <- NULL
  out
}

mid_len = function(mid, sep = ",") {
  mat = strsplit(mid, split = sep, fixed = TRUE)
  mat = do.call(rbind, mat)
  suppressWarnings(storage.mode(mat) <- "numeric")

  out = as.data.frame(cbind(mid, rowSums(!is.na(mat))))
  colnames(out) = c("mid", "size")
  out = out %>% dplyr::mutate(size = as.numeric(size))
  out
}

mid_group = function(mid, sep = ",") {
  mat = strsplit(mid, split = sep, fixed = TRUE)
  mat = do.call(rbind, mat)
  suppressWarnings(storage.mode(mat) <- "numeric")
  mid = mid[order(rowSums(is.na(mat)), decreasing = TRUE)]
  mat = mat[order(rowSums(is.na(mat)), decreasing = TRUE), , drop = FALSE]
  if (nrow(mat) == 1) {
    result = as.data.frame(cbind(mid, mid))
    colnames(result) = c("c", "p")
    return(result)
  }
  mask = as.matrix(is.na(mat))

  NA_flag = (mask %*% t(mask) - rowSums(mask)) == 0
  ones_flag = Longcellsrc::matrix_xor(mat)
  result = Matrix::Matrix(NA_flag & ones_flag, sparse = TRUE)
  result = as.data.frame(Matrix::summary(as(result, "generalMatrix")))
  result = result %>% dplyr::filter(i != j) %>% dplyr::select(-x)

  if (length(result) > 0) {
    result = result[, c("j", "i")]
    colnames(result) = c("c", "p")
    result = result %>% dplyr::mutate(c = mid[c], p = mid[p])
    orphan = setdiff(mid, result$c)
    if (length(orphan) > 0) {
      orphan = as.data.frame(cbind(orphan, NA))
      colnames(orphan) = c("c", "p")
      result = rbind(result, orphan)
    }
  } else {
    result = as.data.frame(cbind(mid, NA))
    colnames(result) = c("c", "p")
  }
  result
}

disagree_sites = function(from, to) {
  if (length(from) != length(to)) {
    stop("There should be a one to one correspondence!")
  }
  from_table = do.call(rbind, strsplit(from, split = ","))
  to_table = do.call(rbind, strsplit(to, split = ","))
  from_table <- suppressWarnings(matrix(as.numeric(from_table), ncol = ncol(from_table)))
  to_table <- suppressWarnings(matrix(as.numeric(to_table), ncol = ncol(to_table)))
  disagree = xor(from_table, to_table) & (from_table == 1)

  filter = as.data.frame(cbind(colSums(disagree, na.rm = TRUE),
                               colSums(from_table, na.rm = TRUE),
                               colSums(to_table, na.rm = TRUE)))
  colnames(filter) = c("disagree", "wrong", "correct")
  filter
}

mid_correct_input = function(cells, cluster, gene_isoform) {
  if (nrow(gene_isoform) != length(cells)) {
    stop("The size of isoforms and cells don't match!")
  }
  if (length(cluster) != length(cells)) {
    stop("The size of clusters and cells don't match!")
  }
  mid = apply(gene_isoform[, 2:(ncol(gene_isoform) - 1), drop = FALSE], 1, function(x) paste(x, collapse = ","))
  out = as.data.frame(cbind(cells, cluster, gene_isoform$start, mid, gene_isoform$end))
  colnames(out) = c("cell", "cluster", "start", "mid", "end")
  out
}

mid_coexist_fast <- function(data) {
  total <- names(sort(table(data$mid), decreasing = TRUE))
  len <- mid_len(total)
  parent <- mid_group(total)

  DT <- data.table::as.data.table(data)
  LEN <- data.table::as.data.table(len)
  data.table::setnames(LEN, c("mid", "size"))
  P <- data.table::as.data.table(parent)
  data.table::setnames(P, c("c", "p"))
  if (nrow(P)) {
    data.table::setkey(P, c, p)
  }

  CNT <- DT[, .N, by = .(cell, cluster, mid)]
  CNT <- LEN[CNT, on = "mid"]
  CNT[is.na(size), size := -Inf]
  data.table::setorder(CNT, cell, cluster, -N, -size)

  TOP2 <- CNT[, head(.SD, 2L), by = .(cell, cluster)]
  NG <- TOP2[, .N, by = .(cell, cluster)]
  groups1 <- TOP2[NG[N == 1L], on = .(cell, cluster)]
  groups2 <- TOP2[NG[N == 2L], on = .(cell, cluster)]

  res_one <- groups1[, .(from = mid, to = mid, count = N)]
  cons_one <- groups1[, .(cell, cluster, concensus = mid)]

  groups2[, other_mid := mid[.N:1L], by = .(cell, cluster)]
  if (nrow(P)) {
    groups2[, is_child_of_other := P[.SD, on = .(c = mid, p = other_mid), .N, by = .EACHI]$N > 0]
  } else {
    groups2[, is_child_of_other := FALSE]
  }

  CONS2 <- groups2[, {
    hd <- sum(is_child_of_other) == 1L
    mm <- if (hd) other_mid[is_child_of_other][1L] else mid[1L]
    sc <- sum(N)
    .(has_dir = hd, mode_mid = mm, sum_count = sc)
  }, by = .(cell, cluster)]

  res_two_collapse <- CONS2[has_dir == TRUE, .(from = mode_mid, to = mode_mid, count = sum_count)]
  res_two_map <- groups2[CONS2[has_dir == FALSE], on = .(cell, cluster)][, .(from = mid, to = mode_mid, count = N)]

  coexist <- data.table::rbindlist(list(res_one, res_two_collapse, res_two_map), use.names = TRUE, fill = TRUE)
  coexist <- coexist[, .(count = sum(as.numeric(count))), by = .(from, to)]

  concensus <- data.table::rbindlist(
    list(cons_one, CONS2[, .(cell, cluster, concensus = mode_mid)]),
    use.names = TRUE, fill = TRUE
  )

  list(concensus, coexist)
}

isoform_corres = function(coexist_matrix) {
  isoform_coexist_filter <- cbind(diag(coexist_matrix), rowSums(coexist_matrix) - diag(coexist_matrix))
  isoform_coexist_filter <- isoform_coexist_filter[order(isoform_coexist_filter[, 1], decreasing = TRUE), ]
  isoform_coexist_filter <- names(which(isoform_coexist_filter[, 1] > 2 * isoform_coexist_filter[, 2]))

  corres = as.data.frame(coexist_matrix[, isoform_coexist_filter, drop = FALSE])
  corres = sapply(seq_len(nrow(corres)), function(i) {
    x = unlist(corres[i, ])
    if (sum(x) == 0) {
      return(NA)
    }
    isoform_coexist_filter[which(x == max(x))[1]]
  })

  corres = as.data.frame(cbind(rownames(coexist_matrix), corres))
  colnames(corres) = c("from", "to")
  rownames(corres) = corres$from

  disagree = corres %>% dplyr::filter(from != to, !is.na(to))
  if (nrow(disagree) == 0) {
    return(corres)
  }
  ds = disagree_sites(disagree$from, disagree$to)
  id = which(ds$wrong > ds$correct & ds$disagree * 3 > ds$wrong)
  if (length(id) == 0) {
    return(corres)
  }

  correct_iso = na.omit(unique(corres$to))
  to_table = do.call(rbind, strsplit(correct_iso, split = ","))
  to_table <- suppressWarnings(matrix(as.numeric(to_table), ncol = ncol(to_table)))
  to_table[is.na(to_table)] = 0
  if (length(id) == 1) {
    correct_iso = correct_iso[to_table[, id] == 0]
  } else if (nrow(to_table) == 1) {
    correct_iso = correct_iso[sum(to_table[, id]) == 0]
  } else {
    correct_iso = correct_iso[rowSums(to_table[, id]) == 0]
  }

  corres_new = as.data.frame(coexist_matrix[, correct_iso, drop = FALSE])
  corres_new = sapply(seq_len(nrow(corres)), function(i) {
    x = unlist(corres_new[i, ])
    if (sum(x) == 0) {
      return(NA)
    }
    correct_iso[which(x == max(x))[1]]
  })

  corres_new = as.data.frame(cbind(rownames(coexist_matrix), corres_new))
  colnames(corres_new) = c("from", "to")
  rownames(corres_new) = corres_new$from
  corres_new
}

cluster_isoform_correct <- function(start, mid, end, concensus, polyA, preserve_mid) {
  if (length(start) != length(mid) || length(start) != length(end)) {
    stop("The size of isoforms representation don't match!")
  }
  if (length(start) != length(polyA)) {
    stop("The size of isoforms and polyA don't match!")
  }

  cluster_size = length(start)
  mode_isoform = NA
  mode_start = NA
  mode_end = NA
  mode_polyA = NA

  if (length(preserve_mid) > 0 && nrow(preserve_mid) > 0) {
    mid_corres = na.omit(preserve_mid[concensus, "to"])
    if (length(mid_corres) > 0) {
      mode_isoform = mid_corres
      mode_preserve = which(mid == mode_isoform)
      if (length(mode_preserve) == 0) {
        mode_start = -1
        mode_end = -1
      } else {
        mode_start = min(start[mode_preserve])
        mode_end = max(end[mode_preserve])
      }
      polyA = as.logical(polyA)
      mode_polyA = mean(as.numeric(polyA[mid == concensus]))
    }
  }

  c(mode_start, mode_isoform, mode_end, cluster_size, mode_polyA)
}

isoform_correct <- function(gene_isoform, preserve_mid) {
  gene_isoform_adjust = gene_isoform %>%
    dplyr::group_by(cell, cluster) %>%
    dplyr::summarise(adjust = list(cluster_isoform_correct(start, mid, end, unique(concensus), polyA, preserve_mid)),
                     .groups = "drop")
  adjust = as.data.frame(do.call(rbind, gene_isoform_adjust$adjust))
  colnames(adjust) = c("start", "mid", "end", "size", "polyA")
  adjust = adjust %>% dplyr::mutate(dplyr::across(c("start", "end", "size", "polyA"), as.numeric))

  out = as.data.frame(cbind(gene_isoform_adjust[, c("cell")], adjust))
  out[!is.na(out$start), ]
}

cells_mid_correct <- function(cells, cluster, gene_isoform, polyA) {
  data = mid_correct_input(cells, cluster, gene_isoform)
  out = mid_coexist_fast(data)
  concensus = as.data.frame(out[[1]])
  coexist = as.data.frame(out[[2]])

  if (nrow(coexist) == 1) {
    corres = as.data.frame(coexist %>% dplyr::select(-count))
    rownames(corres) = corres$from
  } else {
    coexist_matrix = long2square(coexist, "from", "to", "count", symmetric = FALSE)
    corres = isoform_corres(coexist_matrix)
  }

  data = dplyr::left_join(data, concensus, by = c("cell", "cluster"))
  data$polyA = polyA
  isoform_correct(data, corres)
}

cells_nomid_correct <- function(cells, cluster, gene_isoform, polyA) {
  data = as.data.frame(cbind(cells, cluster, gene_isoform, polyA))
  colnames(data) = c("cell", "cluster", "start", "end", "polyA")
  data %>%
    dplyr::group_by(cell, cluster) %>%
    dplyr::summarise(
      start = min(start),
      end = max(end),
      size = dplyr::n(),
      polyA = mean(polyA),
      .groups = "drop"
    )
}

cells_build_isoform_dt <- function(data, sites = NULL, flank = 5L, sep = ",", split = "|") {
  dt <- data.table::as.data.table(data)

  if (is.null(sites)) {
    s_num <- suppressWarnings(as.numeric(dt$start))
    e_num <- suppressWarnings(as.numeric(dt$end))
    dt[, isoform := gsub(" ", "", paste(s_num, e_num, sep = sep), fixed = TRUE)]
    return(dt[])
  }

  sites_num <- suppressWarnings(as.numeric(sites))
  if (anyNA(sites_num)) {
    stop("`sites` must be coercible to numeric.")
  }

  s_num <- suppressWarnings(as.numeric(dt$start))
  e_num <- suppressWarnings(as.numeric(dt$end))

  mid_chr <- as.character(dt$mid)
  mid_chr[is.na(mid_chr)] <- ""
  mids_list <- strsplit(mid_chr, ",", fixed = TRUE)
  mids_num_list <- lapply(mids_list, function(x) if (length(x)) suppressWarnings(as.numeric(x)) else numeric(0L))
  mids_log_list <- lapply(mids_num_list, function(x) suppressWarnings(as.logical(x)))

  if (any(vapply(mids_log_list, length, 1L) != length(sites_num))) {
    stop("The size of splicing sites and binary indicator don't match!")
  }

  masks <- lapply(mids_log_list, function(m) !is.na(m))
  first_idx <- vapply(masks, function(m) if (any(m)) which(m)[1L] else NA_integer_, 1L)
  last_idx <- vapply(masks, function(m) if (any(m)) tail(which(m), 1L) else NA_integer_, 1L)

  repl_s <- !is.na(s_num) & s_num == -1 & !is.na(first_idx)
  repl_e <- !is.na(e_num) & e_num == -1 & !is.na(last_idx)
  if (any(repl_s)) {
    s_num[repl_s] <- sites_num[first_idx[repl_s]] - flank
  }
  if (any(repl_e)) {
    e_num[repl_e] <- sites_num[last_idx[repl_e]] + flank
  }

  dt[, isoform := vapply(seq_len(.N), function(i) {
    si <- s_num[i]
    if (is.na(si)) {
      return(NA_character_)
    }
    ei <- e_num[i]
    m <- mids_log_list[[i]]
    if (length(m)) {
      m[is.na(m)] <- FALSE
    }
    seg <- format(c(si, sites_num[if (length(m)) m else FALSE], ei), scientific = FALSE)
    k <- length(seg) %/% 2L
    if (k == 0L) {
      paste(si, ei, sep = sep)
    } else {
      left <- seg[2L * seq_len(k) - 1L]
      right <- seg[2L * seq_len(k)]
      paste(paste(left, right, sep = sep), collapse = split)
    }
  }, FUN.VALUE = character(1L))]

  dt[, isoform := gsub(" ", "", isoform, fixed = TRUE)]
  dt[]
}

cells_isoform_correct <- function(cells, cluster, gene_isoform, polyA) {
  if (ncol(gene_isoform) > 2) {
    splice_sites = colnames(gene_isoform)[2:(ncol(gene_isoform) - 1)]
    data = cells_mid_correct(cells, cluster, gene_isoform, polyA)
    if (nrow(data) == 0) {
      return(NULL)
    }
    data = cells_build_isoform_dt(data, sites = splice_sites, flank = 5L, sep = ",", split = "|")
  } else {
    data = cells_nomid_correct(cells, cluster, gene_isoform, polyA)
    if (nrow(data) == 0) {
      return(NULL)
    }
    data = cells_build_isoform_dt(data, flank = 5L, sep = ",", split = "|")
    data$mid = "null"
  }

  data = as.data.frame(data)
  data <- na.omit(data) %>% dplyr::select(-c(start, end))
  data
}
