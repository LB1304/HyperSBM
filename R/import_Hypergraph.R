import_hypergraph <- function(file_name, method = c("full", "edge")) {
  if (!file.exists(file_name)) {
    stop("The specified file name does not exist!")
  }
  if (length(method) == 2) {
    method <- method[1]
  }
  
  hg <- data.frame(h_edges = unique(readLines(file_name)))
  hg$sizes <- stringr::str_count(hg$h_edges, ",") + 1
  hg <- hg[order(hg$sizes, decreasing = FALSE), ]
  E <- nrow(hg)
  M <- max(hg$sizes)
  
  List_of_H.Edges <- vector(mode = "list", length = M-1)
  names(List_of_H.Edges) <- paste0(rep("m", M-1), 2:M)
  nodes <- c()
  for (m in unique(hg$sizes)) {
    ind_h_edges_m <- which(hg$sizes == m)
    List_of_H.Edges[[m-1]] <- stringr::str_split_fixed(string = hg$h_edges[ind_h_edges_m], pattern = ",", n = m)
    List_of_H.Edges[[m-1]] <- apply(List_of_H.Edges[[m-1]], 2, as.integer)
    nodes <- sort(unique(c(nodes, unique(as.integer(List_of_H.Edges[[m-1]])))))
  }
  
  # Remove missing nodes (if any)
  while (max(nodes) != length(nodes)) {
    missing_node <- setdiff(1:max(nodes), nodes)[1]
    nodes <- c()
    for (m in unique(hg$sizes)) {
      List_of_H.Edges[[m-1]][List_of_H.Edges[[m-1]] > missing_node] <- List_of_H.Edges[[m-1]][List_of_H.Edges[[m-1]] > missing_node] - 1
      nodes <- sort(unique(c(nodes, unique(as.integer(List_of_H.Edges[[m-1]])))))
    }
  }
  n <- length(nodes)
  
  output <- list(List_of_H.Edges = List_of_H.Edges, Num_nodes = n, Num_H.Edges = E, Max_size = M)
  
  if (method == "full") {
    try(L_full <- Laplacian_full(n = n, E = E, hg = hg, List_of_H.Edges = List_of_H.Edges), silent = TRUE)
    if (exists("L_full")) {
      output <- c(output, L_full)
      output$method <- "full"
    } else {
      L_edge <- Laplacian_edge(n = n, E = nrow(List_of_H.Edges[[1]]), hg = hg, List_of_H.Edges = List_of_H.Edges)
      output <- c(output, L_edge)
      output$method <- "edge"
    }
  } else if (method == "edge") {
    L_edge <- Laplacian_edge(n = n, E = nrow(List_of_H.Edges[[1]]), hg = hg, List_of_H.Edges = List_of_H.Edges)
    output <- c(output, L_edge)
    output$method <- "edge"
  }
  
  class(output) <- "HyperGraph"
  return(output)
}




Laplacian_full <- function(n, E, hg, List_of_H.Edges) {
  H <- matrix(0, nrow = n, ncol = E)
  i <- 1
  for (m in unique(hg$sizes)) {
    mat <- List_of_H.Edges[[m-1]]
    for (e in 1:nrow(mat)) {
      H[mat[e, ], i] <- 1
      i <- i + 1
    }
  }
  D <- rowSums(H)^(-1/2)
  Delta <- colSums(H)^(-1)  #^(-1/2)
  HDelta <- t(t(H) * Delta)
  A <- HDelta %*% t(H)
  L <- diag(n) - (t(D * A) * D)
  
  return(list(Incidence_matrix = H, Nodes_degree = D^(-2), H.Edges_size = Delta^(-2), Laplacian = L))
}


Laplacian_edge <- function(n, E, hg, List_of_H.Edges) {
  H <- matrix(0, nrow = n, ncol = E)
  i <- 1
  mat <- List_of_H.Edges[[1]]
  for (e in 1:nrow(mat)) {
    H[mat[e, ], i] <- 1
    i <- i + 1
  }
  D <- rowSums(H)^(-1/2)
  Delta <- colSums(H)^(-1/2)
  DHDelta <- t(t(D * H) * Delta)
  L <- diag(n) - DHDelta %*% t(DHDelta)
  
  return(list(Incidence_matrix = H, Nodes_degree = D^(-2), H.Edges_size = Delta^(-2), Laplacian = L))
}




