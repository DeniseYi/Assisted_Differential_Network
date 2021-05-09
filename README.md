# Assisted_Differential_Network

Algorithm
-------
Assisted differential network analysis

Maintainer
-------
Huangdi Yi (<huangdi.yi@yale.edu>)


Publication
-------
Yi H, Ma S (2021). Assisted Differential Network Analysis for Gene Expression Data. Manuscript.


Usage
-------
1. ```assisted_differential_network.R```: main functions used in estimation of the proposed assisted differential network analysis.

2. ```networks.R```: functions to generate different network structures and scenarios mentioned in the manuscript.

3. ```main.R```: an example under one simulation setting (p=50, q=50 n=200), containing
   * Generating the true gene networks, regulator networks
   * Generating the corresponding difference matrices
   * A sample simulation setting

Output
------
Two lists containing the estimation for each replicate using the assisted differential network analysis:
  *    ```res_cov_sssvd```: estimation using the covariance matrices as input
  *    ```res_pre_sssvd```: estimation using the precision matrices as input

Four lists containing the estimation for each replicate under using the Alt.1 SVD method: ```res_covX_svd```, ```res_covY_svd```, ```res_preX_svd```, and ```res_preY_svd```.

Four lists containing the estimation for each replicate under using the Alt.2 IRLBA method: ```res_covX_ssvd```, ```res_covY_ssvd```, ```res_preX_ssvd```, and ```res_preY_ssvd```.

Four lists containing the estimation for each replicate under using the Alt.3 bi-clustering SSVD method: ```res_covX_myssvd```, ```res_covY_myssvd```, ```res_preX_myssvd```, and ```res_preY_myssvd```.
