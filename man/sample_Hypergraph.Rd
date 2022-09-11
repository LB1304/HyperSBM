\name{sample_Hypergraph}
\alias{sample_Hypergraph}
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
  \item{Q}{Number of latent groups}
  \item{B}{Estimate of the probabilities that m unordered nodes with a specific latent configuration are connected into a hyperedge}
  \item{pi}{Estimate of the priori probability vector}
  \item{tau}{Matrix of the variational parameters}
  \item{J}{Value of the evidence lower bound J at convergence}
  \item{J_vec}{Vector of the evidence lower bound at every step of the algorithm}
  \item{it}{Number of iterations of the algorithm}
  \item{Z}{Vector containing the estimated latent group for each node}
  \item{LogLik}{Value of the log-likelihood at convergence}
  \item{n_par}{Total number of estimated parameters}
  \item{ICL}{Value of the ICL}
  \item{Y}{List containing the index of the existing hyperedges}
  \item{tau_init}{Matrix of the initial values of the variational parameters}
  \item{call}{Command used to call the function}
}
\examples{
  require(HyperSBM)
  ## Import the co-authorship dataset
  HG <- import_hypergraph(file_name = "HG_coauth.txt", method = "full")
  
  ## Estimation of the HSBM with 2 latent groups, random initialization and full model (not parallelized version)
  est <- HSBM(Hypergraph = HG, Q = 2, tol = 1e-6, maxit_VEM = 100, maxit_FP = 100)
  ## Estimation of the HSBM with 2 latent groups, random initialization and full model (parallelized version with 30 cores)
  est <- HSBM(Hypergraph = HG, Q = 2, tol = 1e-6, maxit_VEM = 100, maxit_FP = 100, n_threads = 30)
  ## Estimation of the HSBM with 2 latent groups, random initialization and affiliation model
  est <- HSBM(Hypergraph = HG, Q = 2, start = 0, model = 2, tol = 1e-6, maxit_VEM = 100, maxit_FP = 100)
  ## Estimation of the HSBM with 2 latent groups, soft spectral clustering initialization and full model
  est <- HSBM(Hypergraph = HG, Q = 2, start = 2, model = 0, tol = 1e-6, maxit_VEM = 100, maxit_FP = 100)
}
\author{Luca Brusa, Catherine Matias}
