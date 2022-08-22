VEM_Q1 <- function(n, M, Q, Y) {
  pi <- 1
  B <- vector(mode = "list", length = M-1)
  for (m in 2:M) {
    B[[m-1]] <- length(Y[[m-1]])/choose(n, m)
  }
  
  tau <- matrix(1, nrow = n, ncol = Q)
  Z <- rep(1, n)
  
  J <- 0
  for (m in 2:M) {
    if (B[[m-1]] == 0 || B[[m-1]] == 1) {
      log_B1 <- 0
      log_B2 <- 0
    } else {
      log_B1 <- log(B[[m-1]])
      log_B2 <- log(1 - B[[m-1]])
    }
    aux <- length(Y[[m-1]]) * log_B1 + (choose(n, m) - length(Y[[m-1]])) * log_B2
    J <- J + aux
  }
  LogLik <- J
  n_par <- M-1
  ICL <- LogLik - 1/2 * sum(log(choose(n, 2:M)))
  
  out = list(Q = Q, B = B, pi = pi, tau = tau, J = J, J_vec = J, it = 1, Z = Z, LogLik = LogLik, n_par = n_par, ICL = ICL, Y = Y, tau_init = tau)
  return(out)
}