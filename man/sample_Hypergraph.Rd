\name{sample_Hypergraph}
\alias{sample_Hypergraph}
\title{Draw samples from the HSBM}
\description{
  This function samples a hypergraph from the affiliation HSBM.
}
\usage{
sample_Hypergraph(n, M, Q, pi, alpha, beta, file_name)
}
\arguments{
  \item{n}{Number of nodes}
  \item{M}{Vector of hyperedges sizes}
  \item{Q}{Number of latent groups}
  \item{pi}{Vector of group probabilities}
  \item{alpha}{Vector of connection probabilities of nodes belonging to the same latent group (same length as M).}
  \item{beta}{Vector of connection probabilities of nodes belonging to at least two different latent groups (same length as M).}
  \item{file_name}{A string with the name of the hypergraph}
}
\details{
  Hypergraphs are sampled under the Aff-m sub-model.
  For each hyperedge size M[i]:
  - alpha[m] represents the probability that M[i] nodes belonging to the same latent block are connected into an hyperedge;
  - beta[m] represents the probability that M[i] nodes belonging to two or more different latent blocks are connected into an hyperedge.
  
  This function has no output but saves two files to the current directory: the "file_name.txt" file, containing the hyperedges, and the "file_name_info.RData" file, containing the vector with the true latent group of each node.
  In the text file, each line represents a different hyperedge, and the nodes in a hyperedge are separated by a comma (e.g. 2, 6, 34).
  
  The hypergraph is drawn from a simplified affiliation model.
}
\examples{
  require(HyperSBM)
  sample_Hypergraph(n = 50, M = c(2, 3), Q = 2, pi = c(0.4, 0.6), alpha = c(0.7, 0.7/50), beta = c(0.3, 0.3/50), file_name = "HG")
}
\author{Luca Brusa, Catherine Matias}
