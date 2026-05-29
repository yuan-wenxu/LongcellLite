frNN_dis <- function(dis, n, eps = 10, if_direct = FALSE) {
  dis <- as.matrix(dis)
  if (!if_direct) {
    dis <- rbind(dis, dis[, c(2, 1, 3)])
  }
  dis <- as.data.frame(dis)
  colnames(dis) <- c("node1", "node2", "dist")
  dis <- dis[order(dis$node1), ]

  out_frNN <- dbscan::frNN(as.dist(1), eps = eps)
  out_frNN$dist <- lapply(seq_len(n), function(i) dis[dis$node1 == i, "dist"])
  out_frNN$id <- lapply(seq_len(n), function(i) dis[dis$node1 == i, "node2"])
  out_frNN$sort <- FALSE
  dbscan::frNN(out_frNN, eps = eps)
}

isoform_dis_cluster <- function(isoforms, thresh = 20, eps = 20, split = "|", sep = ",") {
  isoforms_count = table(isoforms)
  if (length(isoforms_count) == 1) {
    return(rep(1, length(isoforms)))
  }

  dis = Longcellsrc::isos_dis(names(isoforms_count), thresh = thresh, split = split, sep = sep)
  if (nrow(dis) == 0) {
    cluster = seq_along(isoforms_count)
  } else {
    dis_frNN = frNN_dis(dis, length(isoforms_count), eps = eps)
    cluster = dbscan::dbscan(dis_frNN, minPts = 1, weights = isoforms_count)$cluster
  }
  names(cluster) <- names(isoforms_count)
  cluster[isoforms]
}

size_filter_error <- function(size, ratio = 0.1) {
  size <- as.numeric(size)
  if (any(is.na(size))) {
    stop("size must be numeric!")
  }
  n <- length(size)
  if (n <= 1L) {
    return(1)
  }
  if (max(size) == 1) {
    return(rep(1, n))
  }

  ord <- order(size)
  w_sorted <- Longcellsrc::size_filter_cpp(size[ord], ratio)
  w <- numeric(n)
  w[ord] <- w_sorted
  w
}

isoform_size_filter <- function(isoforms, size, ratio = 0.1, ..., thresh = 10, eps = 10) {
  n <- length(isoforms)
  if (n == 0L) {
    return(integer())
  }

  len <- Longcellsrc::isos_len_cpp(isoforms)
  cluster <- isoform_dis_cluster(isoforms, thresh, eps)

  dt <- data.frame(
    cluster = cluster,
    size = as.numeric(size),
    len = as.numeric(len),
    ord = seq_len(n),
    stringsAsFactors = FALSE
  )

  cw <- do.call(rbind, lapply(split(dt, dt$cluster), function(x) {
    data.frame(
      cluster = x$cluster[1],
      weight = sum(size_filter_error(x$size, ratio)),
      stringsAsFactors = FALSE
    )
  }))

  dt <- dt[order(dt$cluster, -dt$size, -dt$len), , drop = FALSE]
  dt$rank <- ave(dt$ord, dt$cluster, FUN = seq_along)
  dt <- merge(dt, cw, by = "cluster", all.x = TRUE, sort = FALSE)
  dt$count <- as.integer(dt$rank <= dt$weight)
  dt <- dt[order(dt$ord), , drop = FALSE]
  dt$count
}

cells_isoforms_size_filter <- function(cell_isoform_table, ratio = 0.1, ..., thresh = 10, eps = 10) {
  dt <- as.data.frame(cell_isoform_table, stringsAsFactors = FALSE)
  group_id <- interaction(dt$cell, dt$mid, drop = TRUE, lex.order = TRUE)
  split_groups <- split(seq_len(nrow(dt)), group_id)

  dt$weight <- 0L
  for (idx in split_groups) {
    dt$weight[idx] <- isoform_size_filter(
      dt$isoform[idx],
      dt$size[idx],
      ratio,
      ...,
      thresh = thresh,
      eps = eps
    )
  }

  out <- do.call(rbind, lapply(split(dt, interaction(dt$cell, dt$isoform, drop = TRUE, lex.order = TRUE)), function(x) {
    sw <- sum(x$weight)
    data.frame(
      cell = x$cell[1],
      isoform = x$isoform[1],
      size = sum(x$size),
      cluster = nrow(x),
      count = sw,
      polyA = if (sw > 0) sum(x$polyA * x$weight) / sw else NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  rownames(out) <- NULL
  out
}
