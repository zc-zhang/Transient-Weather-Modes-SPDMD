% SCALE weather simulation ('datahandle0805.m')
% scaledata240805b.mat (data)  % original SCALE simulation data

% Notations:
% horizontal velocity component in the y-direction (data.V)
% vertical velocity component in the z-direction (data.W) 
% 2D flow in y-z plane corresponding to the velocity componment (data.V, data.W) ==> (TmpV,TmpW)


*Vorticity Magnitude*
% Vfull.mat      % or Vfull09.mat   % magnitude of vorticity data in scalar field by 'vorflow.m' with size 3880*121
% Vfull12.mat    % magnitude of vorticity data, only compute the omega_x (vorticity_x) case by 'vorflow241202.m' 

*Velocity Magnitude*
% XScafull.mat         % magnitude of velocity in scalar field (original data)
% comp_scalar_flow.m   % compute/make a gif for magnitude of scalar flow and generate ('XScafull.mat')
                       % option 1: contourf(.) plot
                       % option 2: imagesc(.) plot


% Other Atmospheric factors:
% pressure (data.PRES) and temperature (data.T) 
% relative humidity (data.QV) and precipitation (data.PREC)


% Pdata.mat     % pressure data
% Tdata.mat     % temperature data
% Hdata.mat     % specific humidity or humidity ratio data


% datahandle0805.m   % handle data (i.e., generating the SCALE data. such as 'Pdata.mat', 'Tdata.ma', 'Hdata.mat') by YS
% vorflow.m          % load('scaledata240805b.mat') to get 'Vfull09.mat' data
% vorflow241202.m    % load('scaledata240805b.mat') to get 'Vfull12.mat' data
