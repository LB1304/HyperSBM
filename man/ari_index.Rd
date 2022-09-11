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
  \item{x}{First clustering}
  \item{y}{Second clustering}
}
\details{
  This function outputs the value of the ARI between two clusterings.
}
\examples{
  require(HyperSBM)
  sample_Hypergraph(n = 50, M = 3, pi = c(0.4, 0.6), alpha = 0.7, beta = 0.3, file_name = "HG")
  
  HG <- import_Hypergraph(file_name = "HG.txt", method = "full")
  res <- HSBM(Hypergraph = HG, Q = 2, )
  
  load("HG_info.RData")
  
  ari_index(x = , y = )
}
\author{Luca Brusa, Catherine Matias}
