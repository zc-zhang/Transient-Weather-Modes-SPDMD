# Transient-Weather-Koopman-Mode-Decomposition
Computation of weather data by means of Koopman mode decomposition (KMD) and its variants. Koopman mode decompostion is an appearing data-driven technique to investigate the complex and nonlinear systems.
Koopman operator theory provides the spectral analysis for the nonlinear behaviors by lifting the nonlinear dynamics into an infinite dimensional but linear system, which acts on the space of the observable functions.
In weather or climate modeling, the true dynamics is usually beyond the mathematical model. However, it is always possible to measure the (short-term) time series data (e.g., SCALE weather simulation data, see: https://scale.riken.jp) .
Based on these weather data at hand, we can exploit data-driven computational method, such as dynamic mode decomposition (DMD) and its variants (e.g., sparsity-promoting dynamic mode decomposition, SPDMD) to further study the spectrum distributions and the corresponding spatial and temporal modes.

# Goal
Extract the transient modes in twerms of warm bubble-like pattern by SPDMD.  


# References:
[1] M. Jovanovic, P. Schmid, and J. Nichols, Sparsity-promoting dynamic mode decomposition, Physics of Fluids, 26, 024103 (2014).  
Doi: https://doi.org/10.1063/1.4863670; dmdsp: https://www.ece.umn.edu/users/mihailo/software/dmdsp/  
[2] S. Nishizawa, H. Yashiro, T. Yamaura, A. Adachi, Y. Sachiho, Y. Sato, and H. Tomita, Scalable Computing for Advanced Library and Environment (SCALE) v5. 3.6 software. 
zenodo, 2020. RIKEN, software: https://scale.riken.jp
