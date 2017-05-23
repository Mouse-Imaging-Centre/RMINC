library(RMINC)
adj <- RMINC:::neighbour_list(50,50,50, 6)

system.time(tfce <- RMINC:::graph_tfce_wqu(mincGetVolume("../../../POND/asymmetry201701101/stats/test_t_downsamp.mnc")
                                           , adj
                                           , E = .5, H = 2, nsteps = 100
                                           , weights = rep(1, 50^3)))

system.time(tfce2 <- vertexTFCE(mincGetVolume("../../../POND/asymmetry201701101/stats/test_t_downsamp.mnc")
                                , igraph::graph.adjlist(lapply(adj, `+`, 1))
                                , E = .5, H = 2, nsteps = 100
                                , side = "positive", weights = rep(1, 50^3)))

system.time(wtfce <- vertexTFCE(mincGetVolume("../../../POND/asymmetry201701101/stats/test_t_downsamp.mnc")
                                , igraph::graph.adjlist(lapply(adj, `+`, 1))
                                , E = .5, H = 2, nsteps = 100
                                , side = "positive"
                                , weight = .001))

system.time(btfce <- vertexTFCE(mincGetVolume("../../../POND/asymmetry201701101/stats/test_t_downsamp.mnc")
                                , igraph::graph.adjlist(lapply(adj, `+`, 1))
                                , E = .5, H = 2, nsteps = 100
                                , side = "both"
                                , weight = .001))

new_adj_list1 <- function(graph, mode = c("all", "in", "out", "total"))
{
  if (!igraph::is_igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph:::igraph.match.arg(mode)
  mode <- as.numeric(switch(mode, out = 1, `in` = 2, all = 3, 
                            total = 3))
  on.exit(.Call("R_igraph_finalizer", PACKAGE = "igraph"))
  res <- .Call("R_igraph_get_adjlist", graph, mode, PACKAGE = "igraph")
  vs <- unclass(igraph::V(graph))
  res <- lapply(res, function(x) vs[x + 1])
  if (igraph::is_named(graph)) 
    names(res) <- igraph::V(graph)$name
  res
}

new_adj_list2 <- function(graph, mode = c("all", "in", "out", "total"))
{
  if (!igraph::is_igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph:::igraph.match.arg(mode)
  mode <- as.numeric(switch(mode, out = 1, `in` = 2, all = 3, 
                            total = 3))
  on.exit(.Call("R_igraph_finalizer", PACKAGE = "igraph"))
  res <- .Call("R_igraph_get_adjlist", graph, mode, PACKAGE = "igraph")
  res <- lapply(res, `+`, 1)
  if (igraph::is_named(graph)) 
    names(res) <- igraph::V(graph)$name
  res
}