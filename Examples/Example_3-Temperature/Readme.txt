% Employ Residual Dynamic Mode Decomposition (Residual DMD) into Temperature data ('Tdata.mat').
% parts of codes stem from MColbrook's GitHub:
https://github.com/MColbrook/Residual-Dynamic-Mode-Decomposition/tree/main

% The file `main-ResDMD' totally stems from MColbrook's GitHub (ref. the above link)

% Test Residual DMD for temperature data to check the point and continuous spectrum.

% Test scaledata using combined Residual DMD and spEDMD (*)  

% Test vortex dynamics (i.e., vorticity) of scaledata using integrated ResidualDMD and spEDMD (**)

(*) spEDMD: sparsity-promoting (extended) dynamic mode decomposition (spEDMD), spDMD
% slightly change the spEDMD form  by taking reweighted L1 norm, or min-max concave penalty to perform spDMD or spEDMD
(**) combined residual and sparsity-promoting DMD (crespDMD)
