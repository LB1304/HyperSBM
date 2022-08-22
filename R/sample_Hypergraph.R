sample_Hypergraph <- function(n, M, Q, pi, alpha, beta, file_name) {
  ## Input:
  #   - n: number of nodes;
  #   - M: maximum size of hyper-edges;
  #   - Q: number of latent blocks;
  #   - pi: vector of weights of each latent block;
  #   - B: conditional probabilities of hyper-edges existence given the configuration of latent blocks.
  
  nodes_in_blocks <- sample(x = 1:Q, size = n, replace = TRUE, prob = pi)
  save(nodes_in_blocks, file = paste0(file_name, "_info.RData"))
  
  sink(paste0(file_name, ".txt"))
  for (m in 2:M) {
    H.Edges_m <- RcppAlgos::comboGeneral(v = n, m = m, repetition = FALSE)
    for (i in 1:nrow(H.Edges_m)) {
      proposed_H.Edge <- H.Edges_m[i, ]
      blocks <- unique(nodes_in_blocks[proposed_H.Edge])
      prob <- ifelse(length(blocks) == 1, alpha, beta)
      exist_H.Edge <- rbinom(n = 1, size = 1, prob = prob)
      if (exist_H.Edge == 1) {
        cat(paste(proposed_H.Edge, collapse = ","), append = TRUE, sep = "\n")
      }
    }
  }
  sink()

}


