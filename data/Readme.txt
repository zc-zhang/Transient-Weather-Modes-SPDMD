% SCALE weather simulation ('datahandle0805.m')
% scaledata240805b.mat (data)  % original SCALE simulation data

% Notations:
% horizontal velocity component in the y-direction (data.V)
% vertical velocity component in the z-direction (data.W) 
% 2D flow in y-z plane corresponding to the velocity componment (data.V, data.W) ==> (TmpV,TmpW)

% Other Atmospheric factors:
% pressure (data.PRES) and temperature (data.T) 
% relative humidity (data.QV) and precipitation (data.PREC)


% Pdata.mat     % pressure data
% Tdata.mat     % temperature data
% Hdata.mat     % specific humidity or humidity ratio data
% Vfull.mat     or Vfull09.mat   % magnitude of vorticity data (scalar field) by 'vorflow.m' with size 3880*121
% Vfull12.mat    % magnitude of vorticityu data, only compute the omega_x (vorticity_x) case by 'vorflow241202.m' 


% XScafull.mat     % magnitude of scale data (scalar filed) does not perform vortex dynamics (original data)

% datahandle0805.m   % handle data (i.e., generating the SCALE data. such as 'Pdata.mat', 'Tdata.ma', 'Hdata.mat') by YS
% vorflow.m          % load('scaledata240805b.mat') to get 'Vfull09.mat' data

