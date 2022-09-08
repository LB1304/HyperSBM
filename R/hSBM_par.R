# Hypergraph: an Hypergraph object
# Q: number of latent blocks
# start: initialization (0 = random, 1 = spectral clustering, 2 = fuzzy spectral clustering)
# model: type of model (0 = full model, 1,2 = affiliation models)

hSBM_par <- function (Hypergraph, Q, M_max = NULL, start = 0, model = 0, tol = 1e-6, maxit_VEM = 50, maxit_FP = 40, n_threads = 1, print = TRUE, seed = NULL) {
  M <- Hypergraph$Max_size
  n <- Hypergraph$Num_nodes
  
  if (class(Hypergraph) != "HyperGraph" || missing(Hypergraph)) {
    stop("An HyperGraph object must be provided!")
  }
  if (missing(Q) || Q <= 0) {
    stop("A correct number of latent blocks must be provided!")
  }
  if (!(start %in% c(0, 1, 2, 3))) {
    stop("A correct initialization method must be provided! (0: random, 1: spectral clustering, 2: fuzzy spectral clustering)")
  }
  if (!(model %in% c(0, 1, 2))) {
    stop("A correct model formulation must be provided! (0: full, 1,2: affiliation)")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(M_max)) {
    M_max <- M
  } else if (!(M_max %in% (3:M))) {
    stop("A correct value for M_max must be provided! (An integer value between 3 and M)")
  }
  
  
  # 1. Compute the list of all possible configurations of latent states
  all_Latents <- vector(mode = "list", length = M-1)
  for (m in 2:M) {
    all_Latents[[m-1]] <- RcppAlgos::comboGeneral(v = Q, m = m, repetition = TRUE)
  }
  
  # 2. Compute Y such that Y[[m]][i] = 1 if the hyper-edge (size = m+2) exists (0 otherwise)
  list_H.Edges <- Hypergraph$List_of_H.Edges
  Y <- compute_Y(n = n, M = M, list_edges = list_H.Edges)
  
  # 3. Caso Q = 1
  if (Q == 1) {
    out <- VEM_Q1(n = n, M = M, Q = Q, Y = Y)
    out$call <- match.call()
    class(out) <- "hSBM"
    return(out)
  }
  
  # 4. Compute the initial values (tau) for the VEM algorithm
  L <- Hypergraph$Laplacian
  if (start == 0) {
    tau_init <- matrix(runif(n*Q), nrow = n, ncol = Q)
    tau_init <- tau_init/rowSums(tau_init)
  } else if (start == 1) {
    X <- as.matrix(eigen(L)$vectors[, (n-Q+1):n])
    X <- X/sqrt(apply(X, 1, function(x) norm(x, "2")))
    X[which(is.nan(X))] <- 0
    km <- kmeans(x = X, centers = Q, nstart = 100)
    tau_vec <- km$cluster
    tau_init <- matrix(0, nrow = n, ncol = Q)
    for (i in 1:n) {
      tau_init[i, tau_vec[i]] <- 1
    }
  } else if (start == 2) {
    X <- as.matrix(eigen(L)$vectors[, (n-Q+1):n])
    X <- X/sqrt(apply(X, 1, function(x) norm(x, "2")))
    X[which(is.nan(X))] <- 0
    f_cm <- ppclust::fcm(x = X, centers = Q, nstart = 100)
    tau_init <- f_cm$u
  } else if (start == 3) {
    eigs <- eigen(L)
    eig_values <- abs(eigs$values)
    Qth <- length(eig_values) - Q
    Qth_eig_value <- sort(eig_values, partial = Qth)[Qth]
    ind_eig <- which(eig_values > Qth_eig_value)
    X <- as.matrix(eigs$vectors[, ind_eig])
    X <- X/sqrt(apply(X, 1, function(x) norm(x, "2")))
    X[which(is.nan(X))] <- 0
    km <- kmeans(x = X, centers = Q)
    tau_vec <- km$cluster
    tau_init <- matrix(0, nrow = n, ncol = Q)
    for (i in 1:n) {
      tau_init[i, tau_vec[i]] <- 1
    }
  }
  
  # 5. VEM algorithm
  out <- VEM_par(M = M_max, Q = Q, n = n, tau_old = tau_init, all_latents = all_Latents[1:M_max-1], Y = Y[1:M_max-1], 
                 tol = tol, maxit_VEM = maxit_VEM, maxit_FP = maxit_FP, model = model, start = start, print = print, n_threads = n_threads)
  Z <- apply(out$tau, 1, which.max)
  out$Z <- Z
  
  # 6. Additional M-Step with complete M_max
  B_complete <- compute_B(n = n, M = M, tau = out$tau, Y = Y, all_latents = all_Latents, model = model, n_threads = n_threads)
  out$B <- B_complete
  pi_complete <- compute_pi(tau = out$tau)
  out$pi <- pi_complete
  
  # 7. ICL criterion
  LogLik <- compute_LogLik(n = n, M = M, Q = Q, tau = out$tau, pi = out$pi, Y = Y, B = out$B)
  out$LogLik <- LogLik
  if (model == 0) {
    n_B_par <- choose(Q+(2:M)-1, 2:M)
    ICL_B_par <- 1/2 * sum(n_B_par * log(choose(n, 2:M)))
  } else if (model == 1) {
    n_B_par <- 2
    ICL_B_par <- sum(log(choose(n, 2:M)))
  } else if (model == 2) {
    n_B_par <- 2 * (M-1)
    ICL_B_par <- (M-1) * sum(log(choose(n, 2:M)))
  }
  out$n_par <- (Q-1) + sum(n_B_par)
  ICL <- LogLik - 1/2 * (Q-1) * log(n) - ICL_B_par
  out$ICL <- ICL
  
  # 8. Output
  out$Y <- Y
  out$tau_init <- tau_init
  out$call <- match.call()
  class(out) <- "hSBM"
  return(out)
}
