coreDetect = function(cores) {
  cores = as.numeric(cores)
  max_cores = parallel::detectCores()
  if (cores > max_cores) {
    warning(paste0(
      "Requested cores exceed the maximum available; using ",
      max_cores,
      " cores instead."
    ))
    cores = max_cores
  }
  cores
}

genes_distribute = function(data, cores, gene_col = "gene") {
  if (cores == 1) {
    return(list(data))
  }

  gene_read_num = data %>%
    group_by(across(all_of(gene_col))) %>%
    summarise(count = n()) %>%
    arrange(desc(count))
  colnames(gene_read_num)[colnames(gene_read_num) == gene_col] = "gene"

  unit = c(seq_len(cores), cores:1)
  gene_read_num = gene_read_num %>%
    mutate(group = rep(unit, (nrow(data) %/% length(unit) + 1))[seq_len(nrow(gene_read_num))])

  gene_list = (gene_read_num %>% group_by(group) %>% summarise(gene = list(gene)))$gene
  lapply(gene_list, function(x) data[data[, gene_col] %in% x, ])
}
