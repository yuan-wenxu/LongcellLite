umi_sim_graph <- function(umi, iso = NULL, sim_thresh = 7, iso_thresh = 200, split = "|", sep = ",") {
  if (is.null(iso)) {
    iso = rep("N", length(umi))
  } else if (length(umi) != length(iso)) {
    warning("The length of umi and isoforms don't correspond, isoform information will be ignored.")
    iso = rep("N", length(umi))
  }

  iso_umi_table = as.data.frame(cbind(iso, umi))
  colnames(iso_umi_table) = c("iso", "umi")
  iso_umi_table$id = seq_len(nrow(iso_umi_table))
  iso_umi_table = iso_umi_table %>%
    dplyr::group_by(iso, umi) %>%
    dplyr::summarise(count = dplyr::n(), id = list(id), .groups = "drop")

  umi_ns = Longcellsrc::umi_graph_table(
    iso_umi_table$umi, iso_umi_table$iso, iso_umi_table$count,
    sim_thresh, iso_thresh, split, sep
  )

  umi_ns = as.data.frame(do.call(rbind, umi_ns))
  colnames(umi_ns) <- c("node1", "node2", "ns", "count")
  umi_ns = umi_ns[umi_ns$ns > 0, ]
  umi_ns$node1 = umi_ns$node1 + 1
  umi_ns$node2 = umi_ns$node2 + 1
  umi_ns = umi_ns %>% dplyr::filter(count > 0) %>% dplyr::mutate(weight = ns)

  graph = igraph::graph_from_data_frame(umi_ns, directed = FALSE, vertices = seq_len(nrow(iso_umi_table)))
  graph = igraph::simplify(graph, remove.loops = TRUE, edge.attr.comb = "first")
  graph = igraph::set_vertex_attr(graph, "count", value = iso_umi_table$count)

  umi_corres_size = sapply(iso_umi_table$id, length)
  umi_corres = rep(seq_len(nrow(iso_umi_table)), umi_corres_size)
  umi_corres = as.data.frame(cbind(do.call(c, iso_umi_table$id), umi_corres))
  colnames(umi_corres) = c("id", "pair_id")
  list(graph, umi_corres)
}

louvain_iter_stack <- function(graph, weight = "weight", resolution = 1, alpha = 2, sim_thresh = 6) {
  in_graph_list = igraph::decompose(graph)
  out_graph_list = list()

  while (length(in_graph_list) > 0) {
    temp_graph = in_graph_list[[1]]
    in_graph_list[[1]] = NULL

    if (length(igraph::V(temp_graph)) == 1) {
      out_graph_list = append(out_graph_list, list(temp_graph))
    } else {
      min_cut <- igraph::graph.mincut(temp_graph)
      if (min_cut >= length(igraph::V(temp_graph)) / alpha ||
          min(igraph::edge_attr(temp_graph)[[weight]]) > sim_thresh) {
        out_graph_list = append(out_graph_list, list(temp_graph))
      } else {
        cluster <- igraph::cluster_leiden(temp_graph, resolution = resolution)
        graph_cut <- igraph::delete_edges(temp_graph, igraph::E(temp_graph)[igraph::crossing(cluster, temp_graph)])
        sub_graphs = igraph::decompose(graph_cut)
        if (length(sub_graphs) == 1) {
          out_graph_list = append(out_graph_list, sub_graphs)
        } else {
          in_graph_list = append(in_graph_list, sub_graphs)
        }
      }
    }
  }
  out_graph_list
}

graph_to_cluster <- function(graph_list) {
  out <- lapply(seq_along(graph_list), function(i) {
    cluster = graph_list[[i]]
    cbind(igraph::vertex_attr(cluster)$name, i)
  })
  out <- as.data.frame(do.call(rbind, out))
  colnames(out) <- c("id", "cluster")
  out$id = as.numeric(out$id)
  out$cluster = as.numeric(out$cluster)
  out <- out[order(out[, "id"]), ]
  out[, "cluster"]
}

umi_cluster <- function(umi, iso = NULL, thresh = NULL) {
  if (is.null(thresh)) {
    thresh = round(length(umi) / 2 + 1)
  }
  graph_corres <- umi_sim_graph(umi, iso = iso, sim_thresh = thresh)
  graph = graph_corres[[1]]
  umi_corres = graph_corres[[2]]

  graph_cluster = louvain_iter_stack(graph = graph, alpha = 2, sim_thresh = thresh + 2)
  cluster = graph_to_cluster(graph_cluster)
  cluster = as.data.frame(cbind(seq_along(cluster), cluster))
  colnames(cluster) = c("pair_id", "cluster")

  cluster_expand = dplyr::left_join(umi_corres, cluster, by = "pair_id")
  cluster_expand = cluster_expand[order(cluster_expand$id), ]
  cluster_expand$cluster
}
