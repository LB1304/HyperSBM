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
  \item{n}{An object of class HyperGraph, as imported by function \code{import_Hypergraph}}
  \item{M}{An integer specifying the number of latent groups}
  \item{Q}{Maximum size of hyperedges}
  \item{pi}{Type of starting value (0=random, 1=spectral clustering, 2=soft spectral clustering, 3=absolute spectral clustreing)}
  \item{alpha}{Model formulation (0=full HSBM mode, 1=, 2=)}
  \item{beta}{Tolerance level for convergence}
  \item{file_name}{Maximum number of iterations of the algorithm}
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
