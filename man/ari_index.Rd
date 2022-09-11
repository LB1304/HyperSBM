\name{ari_index}
\alias{ari_index}
\title{Compute Adjusted Rand Index}
\description{
  This function computes Adjusted Rand Index between two clusterings.
}
\usage{
ari_index(x, y)
}
\arguments{
  \item{x}{Vector with the group of each node in the first clustering}
  \item{y}{Vector with the group of each node in the second clustering}
}
\examples{
  require(HyperSBM)
  sample_Hypergraph(n = 50, M = 3, Q = 2, pi = c(0.4, 0.6), alpha = 0.7, beta = 0.3, file_name = "HG")
  
  HG <- import_Hypergraph(file_name = "HG.txt", method = "full")
  load("HG_info.RData")
  
  res <- HSBM(Hypergraph = HG, Q = 2)
  ari_index(x = nodes_in_blocks, y = res$Z)
}
\author{Luca Brusa, Catherine Matias}
