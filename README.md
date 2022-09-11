<h1 align="center">Stochastic Blockmodel for Hypergraphs</h1>
<p align="center"> <span style="font-size: 14px;"><em><strong>Luca Brusa &middot; Catherine Matias</strong></em></span> </p>
<br>

The `HyperSBM` package contains functions to specify and fit Stochastic Blockmodel for hypergraphs. It allows performing model-based clustering when the relationship among nodes are not limited to pairwise interactions (edges in graphs), but include higher-order interactions among three or more nodes (hyperedges).
The main function relies on a variational implementation of the expectation-maximization algorithm to perform model parameters and latent groups estimation.
It also provides tools for generating hypergraphs according to a simplified latent structure (affiliation sub-model). See the package documentation for a more detailed description of all functions.

To install the `HyperSBM` package directly from GitHub:
```r
# install.packages("devtools")
require(devtools)
devtools::install_github("LB1304/HyperSBM")
```

To download the .tar.gz file (for manual installation) use [this link](https://github.com/LB1304/HyperSBM/archive/main.tar.gz).
