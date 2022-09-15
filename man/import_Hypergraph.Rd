\name{import_Hypergraph}
\alias{import_Hypergraph}
\title{Import hypergraphs}
\description{
  This function imports the hyperedges from a text file and constructs the corresponding hypergraph.
}
\usage{
import_hypergraph(file_name, method = c("full", "edge"))
}
\arguments{
  \item{file_name}{A string with name of the file containing the list of hyperedges.}
  \item{method}{Type of Laplacian computed for the hypergraph: "full" option (default) takes into consideration all the hyperedges, "edge" option only the edges (hyperedges of size 2). The latter could be convenient for computational time purposes.}
}
\details{
  The text file containing the hyperedges must be named as "name.txt".
  Each line of the text file represents a different hyperedge. 
  The nodes in a hyperedge are separated by a comma (e.g. 2,6,34).
}
\value{
  \item{List_of_H.Edges}{List containing the hyperedges}
  \item{Num_nodes}{Total number of nodes}
  \item{Num_H.Edges}{Total number of hyperedges}
  \item{Max_size}{Maximum size of the hyperedges}
  \item{Incidence_matrix}{Incidence matrix}
  \item{Nodes_degree}{Vector of the degree of each node}
  \item{H.Edges_size}{Vector of the size of each hyperedge}
  \item{Laplacian}{Laplacian matrix}
  \item{method}{Type of Laplacian computed for the hypergraph}
}
\examples{
  require(HyperSBM)
  ## Import the co-authorship dataset, computing the full Laplacian matrix
  HG <- import_Hypergraph(file_name = "HG_coauth.txt", method = "full")
  
  ## Import the co-authorship dataset, computing the edge-component Laplacian matrix
  HG <- import_Hypergraph(file_name = "HG_coauth.txt", method = "edge")
}
\author{Luca Brusa, Catherine Matias}
