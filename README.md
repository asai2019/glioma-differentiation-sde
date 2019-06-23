# glioma-differentiation-sde
## Exploring the Information Content of Glioma Differentiation using SDEs
***
### Overview
This repository contains Matlab code to analyze glioma differentiation therapy from an information-theoretic perspective, as part of a publication entitled "Exploring the Information Transmission Properties of Noise-induced Dynamics:
Application to Glioma Differentiation", submitted to BMC Bioinformatics. The files in this repository are partitioned into several experiments meant to quantify the information value of various 

### Model Descrptions
This code applies information theoretic analysis to the Glioma Differentiation 
Network models described by [Sun et al. 2016]:

1. Additive Noise (AN) model: Ito SDEs with constant/additive noise terms
    CT doses = {0, 5, 7.5, 10} (ng/ml)
    noise intensities = {0.1, 1, 5, 10} (%)
(2) Chemical Langevin Equation (CLE) model: Ito SDEs with multiplicative noise terms
    CT doses = {0, 5, 7.5, 10} (ng/ml)
    noise levels = {LL, HL, LH, HH} (intrinsic noise, extrinsic noise)*
(3) Chemical Langevin Equation (CLE-) model with cyclin D1 feedback inhibition:
    This is the same as the CLE model, but with parameter "K6a" multiplied by 10 to
    mimic cyclin D1 positive feedback inhibition
    CT doses = {0, 5, 7.5, 10} (ng/ml)
    noise levels = {LL, HL, LH, HH} (intrinsic noise, extrinsic noise)* 
    *[Noise Setting Values for Intrinsic/Extrinsic Noise]*
    Setting | Standard Deviation of Intrinsic/Extrinsic Noise
    ---------------------------------------------------------
       L                0.001
       H                0.1
----------------------------
Simulation conditions:
Number of Signals: 4 CT doses * 4 noise settings = 16 distinct signals
Number of cells: 500
Duration of simulation: 48 hours
Timestep: 0.01 hours
Resolvable timestep: 1 hour

### Experiments
Experiment | Description
--------- | -----------
0 |  Simulate GFAP levels for 500 simulated cells corresponding to 16 distinct signals, for AN, CLE, and CLE- models
1 |  Compare differentiation potential vs entropy for different noise intensities of AN model 
2 |  Compute heatmaps for summary descriptors (AUC, max response, max fold change) for combinations of CT dose and noise
3 |  Compute heatmaps for summary descriptors for CLE- model
4 |  Compute channel (or information transmission) capacity for summary descriptors for AN, CLE, CLE- models
5 |  Compute channel capacity for AN, CLE, CLE- model for raw and fold-transformed datasets 
6 |  Compute channel capacity for CLE model when asymmetric sampling(balanced, greed)
    a. Balanced sampling: sample equally from each of d subintervals in time domain with maximum variance
    b. Greedy sampling: sample d time points with maximum variance from entire time domain 
7 |  Compute channel capacity for CLE model when clustering response dynamics, and consider only terminally differentiated cells
    (A) Terminally differentiated: cell's GFAP value >= 0.8 (differentiation threshold) at end of simulation
    (B) Clustering response dynamics: apply k-means clustering, with k=3, to entire GFAP profiles for all cells, 
        corresponding to each signal
8 | Compute channel capacity for CLE model when applying principal components analysis to identify d uncorrelated components
      identified from (z-standardized) time series data (d = 1,...,10). This experiment identifies the theoretical limits per se of achieving            optimal information transfer by removing as much correlation as possible, reducing the effects of extrinsic noise in       propagating dependent noise throughout dynamics.








