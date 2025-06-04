# Transient-Weather-Koopman-Mode-Decomposition
Computation of weather data by means of Koopman mode decomposition (KMD) and its variants. Koopman mode decompostion is an appearing data-driven technique to investigate the complex and nonlinear systems.
Koopman operator theory provides the spectral analysis for the nonlinear behaviors by lifting the nonlinear dynamics into an infinite dimensional but linear system, which acts on the space of the observable functions.
In weather or climate modeling, the true dynamics is usually beyond the mathematical model. However, it is always possible to measure the (short-term) time series data (e.g., SCALE weather simulation data, see: https://scale.riken.jp) .
Based on these weather data at hand, we can exploit data-driven computational method, such as dynamic mode decomposition (DMD) and its variants (e.g., sparsity-promoting dynamic mode decomposition, SPDMD) to further study the spectrum distributions and the corresponding spatial and temporal modes.

# Goal
Extract the transient modes in twerms of warm bubble-like pattern by suing SPDMD.  

The matrix \(\mathbf{Y}\) is approximated as:

\[
\mathbf{Y} = \begin{bmatrix}
\mathbf{y}_0 & \mathbf{y}_1 & \cdots & \mathbf{y}_{N-1}
\end{bmatrix}
\approx \mathbf{\Phi}_r \cdot \mathrm{diag}(\mathbf{b}_r) \cdot \mathbf{V}_r
\]

where:
- \(\mathbf{\Phi}_r = \begin{bmatrix} \boldsymbol{\phi}_1 & \boldsymbol{\phi}_2 & \cdots & \boldsymbol{\phi}_r \end{bmatrix}\),
- \(\mathrm{diag}(\mathbf{b}_r) = \begin{bmatrix} b_1 & & & \\ & b_2 & & \\ & & \ddots & \\ & & & b_r \end{bmatrix}\),
- \(\mathbf{V}_r = \begin{bmatrix}
1 & \lambda_1 & \cdots & \lambda_1^{N-1} \\
1 & \lambda_2 & \cdots & \lambda_2^{N-1} \\
\vdots & \vdots & \ddots & \vdots \\
1 & \lambda_r & \cdots & \lambda_r^{N-1}
\end{bmatrix}\).
# References:
[1] M. Jovanovic, P. Schmid, and J. Nichols, Sparsity-promoting dynamic mode decomposition, Physics of Fluids, 26, 024103 (2014).  
Doi: https://doi.org/10.1063/1.4863670; dmdsp: https://www.ece.umn.edu/users/mihailo/software/dmdsp/  
[2] S. Nishizawa, H. Yashiro, T. Yamaura, A. Adachi, Y. Sachiho, Y. Sato, and H. Tomita, Scalable Computing for Advanced Library and Environment (SCALE) v5. 3.6 software. 
zenodo, 2020. RIKEN, software: https://scale.riken.jp
