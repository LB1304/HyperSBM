sample_Hypergraph <- function(n, M, Q, pi, alpha, beta, file_name) {
  ## Input:
  #   - n: number of nodes;
  #   - M: vector of hyper-edges sizes;
  #   - Q: number of latent blocks;
  #   - pi: vector of weights of each latent block;
  #   - alpha: vector of the connection probabilities of nodes in the SAME latent block (same length as M);
  #   - beta: vector of the connection probabilities of nodes in at least two different latent blocks (same length as M).
  #     N.B. Considering alpha and beta the Aff-m and Aff models are implemented; not the full model (for the moment).
  
  
  nodes_in_blocks <- sample(x = 1:Q, size = n, replace = TRUE, prob = pi)
  save(nodes_in_blocks, file = paste0(file_name, "_info.RData"))
  
  sink(paste0(file_name, ".txt"))
  for (m in 1:length(M)) {
    H.Edges_m <- RcppAlgos::comboGeneral(v = n, m = M[m], repetition = FALSE)
    prob <- alpha[m] * apply(H.Edges_m, 1, function(x)length(unique(nodes_in_blocks[x])))
    prob[prob > alpha[m]] <- beta[m]
    exist_H.Edges <- rbinom(prob, 1, prob)
    H.Edges_m <- H.Edges_m[-which(exist_H.Edges == 0), ]
    apply(H.Edges_m, 1, function(x) cat(paste(x, collapse = ","), append = TRUE, sep = "\n"))
  }
  sink()
  
}
