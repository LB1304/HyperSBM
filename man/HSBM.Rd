\name{HSBM}
\alias{HSBM}
\title{Estimation of Hypergraph SBMs}
\description{
  This function performs variational inference of hypergraph Stochastic Block Models (HSBMs).
}
\usage{
HSBM(Hypergraph, Q, M_max = NULL, 
     start = 0, model = 0, 
     tol = 1e-6, maxit_VEM = 50, maxit_FP = 50, 
     n_threads = 1, print = TRUE, seed = NULL)
}
\arguments{
  \item{Hypergraph}{An object of class HyperGraph, as imported by function \code{import_Hypergraph}}
  \item{Q}{An integer specifying the number of latent groups}
  \item{M_max}{Maximum size of hyperedges}
  \item{start}{Type of starting value (0=random, 1=spectral clustering, 2=soft spectral clustering, 3=absolute spectral clustreing)}
  \item{model}{Model formulation (0=full HSBM mode, 1=, 2=)}
  \item{tol}{Tolerance level for convergence}
  \item{maxit_VEM}{Maximum number of iterations of the algorithm}
  \item{maxit_FP}{Maximum number of iterations of the fixed point algorithm}
  \item{n_threads}{Number of cores used in the parallel computation. Default is 1 (not parallelized version)}
  \item{print}{Logical, to print information at each step of the algorithm (default) or not}
  \item{seed}{An integer value with the random number generator state}
}
\value{
  \item{}{}
  \item{}{}
  \item{}{}
  \item{}{}
  \item{}{}
}
\details{
  These are example functions which should be largely
  self-explanatory. Their main benefit is to demonstrate how to write a
  function using the Armadillo C++ classes, and to have to such a
  function accessible from R.
}
\references{
  See the documentation for Armadillo, and RcppArmadillo, for more details.
}
\examples{
  require(HyperSBM)
  HG <- 
  est <- 
}
\author{Luca Brusa, Catherine Matias}

