\name{sample_Hypergraph}
\alias{sample_Hypergraph}
\title{Draw samples from the HSBM}
\description{
  This function draws samples (hypergraphs) from the HSBM.
}
\usage{
sample_Hypergraph(n, M, Q, pi, alpha, beta, file_name)
}
\arguments{
  \item{n}{Number of nodes}
  \item{M}{Maximum size of the hyperedges}
  \item{Q}{Number of latent groups}
  \item{pi}{Vector of priori probabilities}
  \item{alpha}{Probability that nodes in the same latent group are connected into a hyperedge}
  \item{beta}{Probability that nodes in different latent groups are connected into a hyperedge}
  \item{file_name}{A string with the name of the hypergraph}
}
\details{
  This function has no output but saves two files to the current directory: the "file_name.txt" file, containing the hyperedges, and the "file_name_info.RData" file, containing the vector with the true latent group of each node.
  In the text file, each line represents a different hyperedge, and the nodes in a hyperedge are separated by a comma (e.g. 2, 6, 34).
  
  The hypergraph is drawn from a simplified affiliation model.
}
\examples{
  require(HyperSBM)
  sample_Hypergraph(n = 50, M = 3, Q = 2, pi = c(0.4, 0.6), alpha = 0.7, beta = 0.3, file_name = "HG")
}
\author{Luca Brusa, Catherine Matias}
