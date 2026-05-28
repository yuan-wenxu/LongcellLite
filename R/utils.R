saveResult = function(data, filename, sep = "\t") {
  write.table(data, file = filename, sep = sep, quote = FALSE, row.names = FALSE)
}

read2bins = function(read, sep = ",", split = "|") {
  if (is.character(read)) {
    sites = Longcellsrc::isoform2sites(read, split = split, sep = sep)
    len = length(sites)
    bins = as.data.frame(cbind(sites[seq(1, len, 2)], sites[seq(2, len, 2)]))
  } else if (is.matrix(read) || is.data.frame(read)) {
    if (ncol(read) != 2) {
      stop("There should be exactly 2 columns for a bin matrix.")
    }
    bins = as.data.frame(read)
  } else {
    stop("Unsupported input type for read2bins().")
  }

  colnames(bins) = c("start", "end")
  bins
}

long2square = function(long, row_names_from, col_names_from, values_from,
                       symmetric = TRUE, na.fill = 0, nodes = NULL) {
  long = long[, c(row_names_from, col_names_from, values_from)]
  if (symmetric) {
    rlong = long[, c(col_names_from, row_names_from, values_from)]
    colnames(rlong) = c(row_names_from, col_names_from, values_from)
    long = rbind(long, rlong)
    long = long[!duplicated(long), ]
  }
  if (is.null(nodes)) {
    nodes = unique(unlist(c(long[, row_names_from], long[, col_names_from])))
  }
  mat = as.data.frame(tidyr::pivot_wider(long, names_from = col_names_from, values_from = values_from))
  rownames(mat) = mat[, row_names_from]
  mat = mat[, -which(colnames(mat) == row_names_from), drop = FALSE]
  diff = setdiff(nodes, colnames(mat))
  if (length(diff) > 0) {
    mat[, diff] = NA
  }
  mat = mat[nodes, nodes, drop = FALSE]
  rownames(mat) = colnames(mat) = nodes
  mat[is.na(mat)] = na.fill
  as.matrix(mat)
}

save10X = function(long, path, i = "gene", j = "cell", value = "count") {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  long = long[, c(i, j, value)]
  x = names(table(long[, i]))
  y = names(table(long[, j]))
  long = long %>% dplyr::mutate(dplyr::across(c(i, j), ~ as.numeric(as.factor(.))))
  long = as.data.frame(long)

  sparse_mat = Matrix::sparseMatrix(
    i = long[, i],
    j = long[, j],
    x = long[, value]
  )

  write.table(x, file = file.path(path, "features.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(y, file = file.path(path, "barcodes.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  Matrix::writeMM(sparse_mat, file = file.path(path, "matrix.mtx"))

  R.utils::gzip(file.path(path, "features.tsv"), overwrite = TRUE)
  R.utils::gzip(file.path(path, "barcodes.tsv"), overwrite = TRUE)
  R.utils::gzip(file.path(path, "matrix.mtx"), overwrite = TRUE)
}

saveIsoMat = function(iso, path, cell_col = "cell", gene_col = "gene",
                      iso_col = "isoform", count_col = "count") {
  dir.create(file.path(path, "gene"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(path, "isoform"), showWarnings = FALSE, recursive = TRUE)

  iso = iso[, c(cell_col, gene_col, iso_col, count_col)]
  colnames(iso) = c("cell", "gene", "isoform", "count")

  gene_long = iso %>% dplyr::group_by(cell, gene) %>%
    dplyr::summarise(count = sum(count), .groups = "drop")
  iso_long = iso %>% dplyr::filter(isoform != "unknown", count > 0) %>%
    dplyr::select(-gene)

  save10X(gene_long, file.path(path, "gene"))
  save10X(iso_long, file.path(path, "isoform"), i = "isoform")
}

init_project = function(work_dir) {
  dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(work_dir, "annotation"), showWarnings = FALSE)
  dir.create(file.path(work_dir, "bam"), showWarnings = FALSE)
  dir.create(file.path(work_dir, "read_isoform"), showWarnings = FALSE)
  dir.create(file.path(work_dir, "out"), showWarnings = FALSE)
  normalizePath(work_dir)
}
