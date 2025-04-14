% plot the spatial modes and temporal 
% plot_STmodes_spdmd.m

%% sort amplitudes for xpol, xsp, xdmd
% command it (if have been done)
% 
 [Norm_xdmd,Index_xdmd] = sort(abs(xdmd),'descend');
 DEv_xdmd = Edmd(Index_xdmd);   %discrete-eigenvalues 
 DMDModes_xdmd = Phi(:,Index_xdmd);

 
 kk=50;
% sort of the large amplitudes rather than rank of Eigvals
 [Norm_xsp,Index_xsp] = sort(abs(answer.xsp(:,kk)),'descend');  %return the value of order/peak ofamplitudes
 DEv_xsp = Edmd(Index_xsp);   %discrete-eigenvalues 
 DMDModes_xsp = Phi(:,Index_xsp);
 

 [Norm_xpol,Index_xpol] = sort(abs(answer.xpol(:,kk)),'descend');
 DEv_xpol = Edmd(Index_xpol);   %discrete-eigenvalues 
 DMDModes_xpol = Phi(:,Index_xpol);


%% tempporal modes 
% Define modes to plot (example: 1, 3, 5, 7, 9, 11)
% kk =50;
rr=answer.Nz(kk);
modesToPlot = 1:2:rr; % Adjust this based on your actual mode range
numModes = length(modesToPlot);

% Number of rows and columns in the subplot grid
numRows = 2;
numCols = 3;

% Define the row and column of the point you want to highlight, in[40,97]
 highlight_row = 26;  % <40 
 highlight_col = 27;   % 97

% Calculate the linear index corresponding to this location
highlight_index = (highlight_col - 1) * 40 + highlight_row;

% Create a figure for the subplots
figure;

for k = 1:numModes
    i = modesToPlot(k);
    
    % Determine subplot position
    subplotRow = ceil(k / numCols);
    subplotCol = mod(k - 1, numCols) + 1;
    
    % Extract and reshape mode_i with amplituide (resp., real, image, angle parts)
    %mode_i = abs(DMDModes_xdmd(:, i)); % DMD modes
    %mode_i = abs(CKModes(:, i));       % Companian Modes
    mode_i = abs(DMDModes_xsp(:, i));  % spDMDModes
    mode_i = mode_i(1:40*97); % Extract the first 40*97 elements
    
    % Plot the spatial structure in the subplot
    subplot(numRows, numCols, (subplotRow - 1) * numCols + subplotCol);
    imagesc(data.y, data.z, reshape(mode_i, [40, 97])'); % Use transpose for correct orientatio  
    colorbar; % Add colorbar for reference
    colormap(brighten(redblueTecplot(21),-0.55));
    axis xy; % Ensure the correct orientation
    xlim([0 2e4]);
    ylim([0 2e4]);
    xlabel("y");
    ylabel("z"); 
    title(['Spatial Mode ', num2str(i)]);
    
    % Highlight the specific point
    hold on;
    % Calculate the actual coordinates for y and z
    highlight_y = data.y(highlight_row);
    highlight_z = data.z(highlight_col);
    plot(highlight_y, highlight_z, 'ko', 'MarkerSize', 9, 'MarkerFaceColor', 'yellow'); % Use a yellow circle
    hold off;
end




%===============================================================
%%%% temporal modes of selsted + the i-the element of recosntraucted Vspdmd
% Define time vector and initialize time dynamics
t = linspace(0, 42, size(V0, 2));  % user-specific time horizion
time_dynamics = zeros(r, length(t));

% Compute the time dynamics for each time step
for iter = 1:length(t)
    %time_dynamics(:, iter) = Norm_xsp .* exp(log(DEv_xsp) * t(iter));
     time_dynamics(:, iter) = answer.xsp(:,kk) .* exp(log(DEv_xsp) * t(iter));
end

% Compute Vspdmd
% Vspdmd = DMDModes_xsp * time_dynamics;
Vspdmd = Phi *diag(answer.xsp(:,kk))*Vand  % Phi =DMDModes_xsp;

%save('Vspdmd',"Vspdmd")


%% Plot Eigenvalues (spDMD vs DMD)
nonzero_indices = find(answer.xsp(:,kk)); % Step 1: Find nonzero amplitudes
[~, sort_order] = sort(nonzero_indices, 'descend'); 
ival = nonzero_indices(sort_order);    % also ralted to `DEv_xsp' (the order of eigvals)
% Spectrum of DT system via SpDMD vs DMD with different Nz
figure;
plot(Edmd, 'bo');
hold on
plot(Edmd(ival), 'r+');
rectangle('Position', [-1 -1 2 2], 'Curvature', 1, ...
'EdgeColor', 'k', 'LineStyle', '--');
xlabel('Real part');
 ylabel('Imaginary part');
 title('Eigenvalues SCALE')
axis (1.2*[-1 1 -1 1]);
axis square;
grid on

% Create a single figure for combined plots
figure;

%highlight_index=1190;
%  highlight_row = 26;
%  highlight_col = 27;
% 
% % Calculate the linear index corresponding to this location
% highlight_index = (highlight_col - 1) * 40 + highlight_row;
% Plot the temporal dynamics for the specific highlight index
plot(t, real(Vspdmd(highlight_index, :)), 'LineWidth', 2, 'DisplayName', ['Highlight Index ', num2str(highlight_index)]);
hold on;

plot(t, real(V0(highlight_index, :)), 'LineWidth', 2, 'DisplayName', ['Original Index ', num2str(highlight_index)]);
hold on;

% Plot the temporal dynamics for selected modes
for i = 1:2:rr
    plot(t, real(time_dynamics(i, :)), 'LineWidth', 2, 'DisplayName', ['Temporal Mode ', num2str(i)]);
end

% Customize the plot
title('Temporal Dynamics Comparison');
xlabel('Time');
ylabel('Evolution');
legend('show'); % Automatically uses the 'DisplayName' property to create legend entries
grid on;
hold off;


%%% Comparison original DMD-temporal

t = linspace(0, 42, size(V0, 2));
time_dynamics_Org = zeros(r, length(t));

% Compute the time dynamics for each time step
for iter = 1:length(t)
    %time_dynamics_Org(:, iter) = Norm_xdmd .* exp(log(DEv_xdmd) * t(iter));
     %  time_dynamics_Org(:, iter) = Norm_xdmd  .* exp(log(DEv_xdmd) * t(iter));
      time_dynamics_Org(:, iter) = xdmd  .* exp(log(DEv_xdmd) * t(iter));
end

% Compute Vspdmd
%Vdmd = DMDModes_xdmd * time_dynamics_Org;
Vdmd = Phi * time_dynamics_Org;
% Create a single figure for combined plots

%% compute Companion DMD
 [CKModes,CKEv,CNorms,VCdmd] = CompanionMatrix_DMD( Vfull );


%option: Residual DMD testing
% Residual DMD
% [G,Kres,Lres,PX,PY] = kernel_ResDMD(V0,V1,'type',"Gaussian");
% [ResW,LAMres,W2] = eig(Kres,'vector');
%  Res = abs(sqrt(real(diag(W2'*Lres*W2)./diag(W2'*W2)-abs(LAMres).^2)));



figure
%scatter(real(LAMres),imag(LAMres),250,R,'.','LineWidth',1);
%hold on
scatter(real(LAMres),imag(LAMres),'.'); hold on
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'-k')
axis equal
axis([-1.15,1.15,-1.15,1.15])
clim([0,1])
 colorbar;
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18)
%title(sprintf('Residuals ($M=%d$)',M),'interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

%exportgraphics(gcf,sprintf('Vlow_res_M%d.pdf',M),'ContentType','vector','BackgroundColor','none')

figure
loglog([0.001,1],[0.001,1],'k','linewidth',2)
hold on
loglog(sqrt(abs(abs(LAMres).^2-1)),R,'b.','markersize',20)
xlabel('$\sqrt{|1-|\lambda|^2|}$','interpreter','latex','fontsize',18)
ylabel('residual','interpreter','latex','fontsize',18)
%title(sprintf('Residuals ($M=%d$)',M),'interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

exportgraphics(gcf,sprintf('Vflow_res2_M%d.pdf',M),'ContentType','vector','BackgroundColor','none')


 %% Pyhsical-informed DMDM (piDMD)
 % Apidmd = piDMD(V0,V1,'diagonaltls');

figure;

% %highlight_index=1190;
%   highlight_row = 24;
%   highlight_col = 19;
% 
% % % Calculate the linear index corresponding to this location
%  highlight_index = (highlight_col - 1) * 40 + highlight_row;

% Plot the temporal dynamics for the specific highlight index
plot(t, real(V0(highlight_index, :)), 'LineWidth', 2, 'DisplayName', ['Original Index ', num2str(highlight_index)]);
hold on;

plot(t, real(VCdmd(highlight_index, :)), 'LineWidth', 2, 'DisplayName', ['CDMD-Recon y_k,i ', num2str(highlight_index)]);
hold on;

plot(t, real(Vdmd(highlight_index, :)), 'LineWidth', 2, 'DisplayName', ['DMD-Recon y_k,i ', num2str(highlight_index)]);
hold on;

plot(t, real(Vspdmd(highlight_index, :)), 'LineWidth', 2, 'DisplayName', ['SPDMD-Recon y_k,i ', num2str(highlight_index)]);
hold on;


% Customize the plot
title('Temporal Dynamics Comparison');
xlabel('Time');
ylabel('Evolution');
legend('show'); % Automatically uses the 'DisplayName' property to create legend entries
grid on;
hold off;



%% Plot Root Mean square error (RMSE)
% Compute RMSE over time for DMD-reconstructed data in specific loacation
rmse_Vdmd = sqrt(mean((real(V0(highlight_index, :) - Vdmd(highlight_index, :))).^2, 1))/ sqrt(length(t));

% Compute RMSE over time for SPDMD-reconstructed data
rmse_Vspdmd = sqrt(mean((real(V0(highlight_index, :) - Vspdmd(highlight_index, :))).^2, 1))/ sqrt(length(t));

% Plot the RMSE over time for DMD and SPDMD
figure;
plot(t, rmse_Vdmd, 'LineWidth', 2, 'DisplayName', 'RMSE DMD');
hold on;
plot(t, rmse_Vspdmd, 'LineWidth', 2, 'DisplayName', 'RMSE SPDMD');
xlabel('Time');
ylabel('RMSE');
title(['Root Mean Square Error ', num2str(highlight_index)]);
legend('show');
grid on;
hold off;



% Case 1: V0
mean_V0 = mean(real(V0), 1);   % Mean over all 3880 trajectories, for each time step
std_V0 = std(real(V0), 0, 1);  % Standard deviation over all 3880 trajectories, for each time step

% Case 2: Vdmd
mean_Vdmd = mean(real(Vdmd), 1);   % Mean for Vdmd over 3880 trajectories
std_Vdmd = std(real(Vdmd), 0, 1);  % Standard deviation for Vdmd

% Case 3: VCdmd
mean_VCdmd = mean(real(VCdmd), 1);  % Mean for VCdmd over 3880 trajectories
std_VCdmd = std(real(VCdmd), 0, 1);  % Standard deviation for VCdmd

% Case 4: Vspdmd
mean_Vspdmd = mean(real(Vspdmd), 1);   % Mean for Vspdmd over 3880 trajectories
std_Vspdmd = std(real(Vspdmd), 0, 1);  % Standard deviation for Vspdmd

% Plot Case 1: V0
figure;
shadedErrorBar(t, mean_V0, std_V0, 'lineprops', '--b', 'patchSaturation', 0.33);
hold on;

% Plot Case 2: Vdmd
shadedErrorBar(t, mean_Vdmd, std_Vdmd, 'lineprops', {'-go', 'MarkerFaceColor', 'g'}, 'patchSaturation', 0.33);

% Plot Case 3: VCdmd
shadedErrorBar(t, mean_VCdmd, std_VCdmd, 'lineprops', '-r', 'patchSaturation', 0.33);

% Plot Case 4: Vspdmd
shadedErrorBar(t, mean_Vspdmd, std_Vspdmd, 'lineprops', {'-mo', 'MarkerFaceColor', 'm'}, 'patchSaturation', 0.33);

% Customizations
xlabel('Time');
ylabel('Amplitude');
title('Comparison of Four Trajectories (V0, Vdmd, VCdmd, Vspdmd) with Shaded Error Bars');
grid on;
legend({'V0', 'Vdmd', 'VCdmd', 'Vspdmd'}, 'Location', 'Best');
hold off;




