# HLloyd
This is the R package for paper: "Exact Clustering in Tensor Block Model: Statistical Optimality and Computational Limit" by Rungang Han, Yuetian Luo, Miaoyan Wang and Anru Zhang (2020).

# Instructions
HLloyd requires the following packages for full functionality: 'rTensor', 'gtools', 'LICORS', 'mclust'. Use the following commands for installation:

library(devtools)  
devtools::install_github("Rungang/HLloyd")

See example.R for illustration.

Please contact hrg.stat@gmail.com with any questions or comments.

# Citation
@article{han2020exact,
  title={Exact Clustering in Tensor Block Model: Statistical Optimality and Computational Limit},
  author={Han, Rungang and Luo, Yuetian and Wang, Miaoyan and Zhang, Anru R},
  journal={arXiv preprint arXiv:2012.09996},
  year={2020}
}

# Reproducibility
The code for obtaining numeric results are collected in the /experiment subfolder.

- All main functions are in HOLloyd.R 

- Simulation codes in Section 6.1:
  - matrix_phase_transition.R and tensor_phase_transition.R: Statistical and Computational Phase Transition.
  - different_delta.R: Impact of Initialization to HLloyd.
  - tensor_recovery.R: Clustering via HLloyd + HSC algorithms.
  - tensor_recovery.R: Tensor estimation via HLloyd + HSC algorithms.
  - bernoulli_model_clustering_property.R: Simulations on Stochastic Tensor Block Models.

- Simultion codes in Section 6.2:
  - TC_methods_compare.R, TC_methods_compare_diffr.R, TC_methods_compare_imbalance.R: Three comparisons with baseline methods.

- Real data in Section 7.1 (flight route network)
  - Data files: flight_route.RData, US_flight_route.RData
  - code files: airline_clustering.R and US_airline_clustering.R
  - Notice: In the simulatin code, the 'delta' parameter in the code relates to the 'gamma' parameter in the paper with relation '-2 delta = gamma' 

- Real data in Section 7.2 (Online Click-through Data)
  - Data file: Click_through.csv
  - code files: click_through.R
  - Notice: In the data file, the four columns represent the user ID, item ID, timepoint (1~24h) and date (1~8 day).
