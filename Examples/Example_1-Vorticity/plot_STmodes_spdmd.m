% plot the spatial modes and temporal 
% plot_spdmd_vorflow
% plot_STmodes_spdmd.m
clc;
% % 
[Norm_xdmd,Index_xdmd] = sort(abs(xdmd),'descend');
DEv_xdmd = Edmd(Index_xdmd);   %discrete-eigenvalues 
xdmd_sorted = xdmd(Index_xdmd);
DMDModes_xdmd = Phi(:,Index_xdmd);

%kk =50; %vorticity
kk=87;  %kk=123;
rr=answer.Nz(kk);

%sort_xsp = answer.xsp(:,kk);  
%sort of the large amplitudes rather than rank of Eigvals
[Norm_xsp,Index_xsp] = sort(abs(answer.xsp(:,kk)),'descend');  %return the value of order/peak of amplitudes
DEv_xsp = Edmd(Index_xsp);   %discrete-eigenvalues 
xsp_sorted = answer.xsp(Index_xsp,kk);
DMDModes_xsp = Phi(:,Index_xsp);


[Norm_xpol,Index_xpol] = sort(abs(answer.xpol(:,kk)),'descend');
DEv_xpol = Edmd(Index_xpol);   %discrete-eigenvalues 
xpol_sorted = answer.xpol(Index_xpol,kk);
DMDModes_xpol = Phi(:,Index_xpol);


%%   Eigenvalues bar by reordered of its amplitude (absolute value of eigenvalues)
% Step 1: Sort the magnitudes of Edmd in descending order
[magnitudes, sorted_indices] = sort(abs(Edmd), 'descend');  % Sort magnitudes in descending order
sorted_Edmd = Edmd(sorted_indices);                         % Reorder Edmd according to sorted magnitudes

% Step 2: Find positions of Edmd(ival) in the sorted list for color overlay
[~, overlay_indices] = ismember(ival, sorted_indices);  % Find positions of Edmd(ival) in sorted order
overlay_indices = overlay_indices(overlay_indices > 0); % Remove any zero entries (if any indices do not match)

% Step 3: Plot all magnitudes with different colors for the Edmd(ival) subset
figure;
hold on;

% Plot all magnitudes in a default color (e.g., light gray)
bar(1:length(magnitudes), magnitudes, 'FaceColor', [0.8, 0.8, 0.8]);

% Overlay bars corresponding to Edmd(ival) in blue
bar(overlay_indices, magnitudes(overlay_indices), 'FaceColor', 'b'); % Blue for Edmd(ival) bars

% Step 4: Add title and labels
title('Magnitude of Eigenvalues with Edmd(ival) Highlighted');
xlabel('Eigenvalue Index (sorted by descending magnitude)');
ylabel('Magnitude');

% Optional legend for clarity
legend({'Other Magnitudes', 'Edmd(ival) Highlighted'}, 'Location', 'Best');

% Add grid and release hold
grid on;
hold off;

%% spatial modes 
% Define the modes to plot and highlight location
modesToPlot = 1:2:rr;
numModes = length(modesToPlot);

highlight_row = 25;    % In [1, 40]
highlight_col = 26;    % In [1, 97]
highlight_index = (highlight_col - 1) * 40 + highlight_row;

% Set up figure and tiled layout
figure;
numRows = 2;
numCols = 3;
ttt = tiledlayout(numRows, numCols, 'Padding', 'compact', 'TileSpacing', 'compact');
%title(t, 'Dominant Spatial Modes', 'FontSize', 22, 'FontWeight', 'bold');

maxTiles = numRows * numCols;
plotCount = min(numModes, maxTiles);

for k = 1:plotCount
    i = modesToPlot(k);
    mode_i = abs(DMDModes_xsp(:, i));
    mode_i = mode_i(1:40 * 97);
    mode_i_reshaped = reshape(mode_i, [40, 97])';

    ax = nexttile;
    imagesc(data.y, data.z, mode_i_reshaped);
    colormap(ax, brighten(redblueTecplot(21), -0.55));
    colorbar;
    
    axis xy;
    xlim([0 2e4]);
    ylim([0 2e4]);
    set(gca, 'FontSize', 17.6);  % Set tick label font size
    xlabel("y", 'FontSize', 17.6);
    ylabel("z", 'FontSize', 17.6);
    title(['Spatial Mode ', num2str(i)], 'FontSize', 17);

    % Highlight specific point
    highlight_y = data.y(highlight_row);
    highlight_z = data.z(highlight_col);
    hold on;
    plot(highlight_y, highlight_z, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'yellow');
    hold off;
end

% Optional: save figure
set(gcf, 'Renderer', 'painters');


%%%% temporal modes of selsted + the i-the element of recosntraucted Vspdmd
% Define time vector and initialize time dynamics


% Compute the time dynamics for each time step
r=rank(Fdmd);
 t = linspace(0, r, size(V0, 2));
time_dynamics = zeros(r, length(t));       % SPDMD reconstruted 
time_dynamics_Org = zeros(r, length(t));   % DMD reconstrauted 

% Compute the time dynamics for each time step
for iter = 1:length(t)
    time_dynamics(:, iter) = xsp_sorted.* exp(log(DEv_xsp) * t(iter));
    time_dynamics_Org(:, iter) = xdmd_sorted  .* exp(log(DEv_xdmd) * t(iter));
end


% Compute Vspdmd
%Vspdmd = DMDModes_xsp * time_dynamics;
Vspdmd=Phi*diag(answer.xsp(:,kk))*Vand;

% cpompute Vdmd
Vdmd = DMDModes_xdmd *time_dynamics_Org;
% Vdmd =Phi*diag(xdmd)*Vand;


%save('Vspdmd',"Vspdmd")


 %% plot the temoral dyanmics
% Define modes to plot (same as spatial)
modesToPlot = 1:2:rr;
numModes = length(modesToPlot);

% Set up tiled layout
figure;
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for k = 1:min(numModes, 8)  % Avoid plotting more than 8
    i = modesToPlot(k);
    nexttile;

    % Plot temporal dynamics for mode i
    plot(t, real(time_dynamics(i, :)), 'LineWidth', 2, 'Color', [0 0.447 0.741]);
    grid on;

    % Customize axes
    xlabel('Time', 'FontSize', 17.6);
    ylabel('Evolution', 'FontSize', 17.6);
    title(['Temporal Mode ', num2str(i)], 'FontSize', 17);
    set(gca, 'FontSize', 13);
end

% Create a figures for combined plots
%figure;

%highlight_index=1190;
%  highlight_row = 22;
%  highlight_col = 13;
% 
% % Calculate the linear index corresponding to this location
% highlight_index = (highlight_col - 1) * 40 + highlight_row;
% Plot the temporal dynamics for the specific highlight index
% plot(t, real(Vspdmd(highlight_index, :)), 'LineWidth', 2, 'DisplayName', ['Highlight Index ', num2str(highlight_index)]);
% hold on;

% plot(t, real(V0(highlight_index, :)), 'LineWidth', 2, 'DisplayName', ['Original Index ', num2str(highlight_index)]);
 %hold on;






%% compute Companion DMD
% [CKModes,CKEv,CNorms,VCdmd] = CompanionMatrix_DMD( Vfull );



 %% Residual DMD
 %[G,K,L,PX,PY,PSI_x,PSI_y,PSI_y2] = kernel_ResDMD(V0,V1,'type',"Gaussian");
%  [G,Kres,Lres,PX,PY] = kernel_ResDMD(V0,V1,'type',"Gaussian");
%  [ResW,LAMres,W2] = eig(Kres,'vector');
% Res = abs(sqrt(real(diag(W2'*Lres*W2)./diag(W2'*W2)-abs(LAMres).^2)));
% 
% figure
% %scatter(real(LAMres),imag(LAMres),250,R,'.','LineWidth',1);
% %hold on
% scatter(real(LAMres),imag(LAMres),'.'); hold on
% plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'-k')
% axis equal
% axis([-1.15,1.15,-1.15,1.15])
% clim([0,1])
%  colorbar;
% xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18)
% ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18)
% %title(sprintf('Residuals ($M=%d$)',M),'interpreter','latex','fontsize',18)
% ax=gca; ax.FontSize=18;
% 
% %exportgraphics(gcf,sprintf('Vlow_res_M%d.pdf',M),'ContentType','vector','BackgroundColor','none')
% 
% figure
% loglog([0.001,1],[0.001,1],'k','linewidth',2)
% hold on
% loglog(sqrt(abs(abs(LAMres).^2-1)),R,'b.','markersize',20)
% xlabel('$\sqrt{|1-|\lambda|^2|}$','interpreter','latex','fontsize',18)
% ylabel('residual','interpreter','latex','fontsize',18)
% %title(sprintf('Residuals ($M=%d$)',M),'interpreter','latex','fontsize',18)
% ax=gca; ax.FontSize=18;

%exportgraphics(gcf,sprintf('Vflow_res2_M%d.pdf',M),'ContentType','vector','BackgroundColor','none')


 %% Pyhsical-informed DMDM (piDMD)
 % Apidmd = piDMD(V0,V1,'diagonaltls');

figure;

% % % Calculate the linear index corresponding to this location
%  highlight_index = (highlight_col - 1) * 40 + highlight_row;

% Plot the temporal dynamics for the specific highlight index
plot(t, abs(V0(highlight_index, :)), 'LineWidth', 2, 'DisplayName', ['Original Index ', num2str(highlight_index)]);
hold on;

% plot(t, abs(VCdmd(highlight_index, :)), 'LineWidth', 2, 'DisplayName', ['CDMD-Recon y_k,i ', num2str(highlight_index)]);
% hold on;

plot(t, abs(Vdmd(highlight_index, :)), 'LineWidth', 2, 'DisplayName', ['DMD-Recon y_k,i ', num2str(highlight_index)]);
 hold on;

plot(t, abs(Vspdmd(highlight_index, :)), 'LineWidth', 2, 'DisplayName', ['SPDMD-Recon y_k,i ', num2str(highlight_index)]);
hold on;


% Customize the plot
title('Temporal Dynamics Comparison');
xlabel('Time');
ylabel('Evolution');
legend('show'); % Automatically uses the 'DisplayName' property to create legend entries
grid on;
hold off;

%% optional: realx the optional
% Case 1: oiginal Humidity ratio data V0
mean_V0 = mean(real(V0), 1);   % Mean over all 3880 trajectories, for each time step
std_V0 = std(real(V0), 0, 1);  % Standard deviation over all 3880 trajectories, for each time step

% Case 3: DMD reconstruction
mean_Vdmd = mean(real(Vdmd), 1);   % Mean for Vdmd over 3880 trajectories
std_Vdmd = std(real(Vdmd), 0, 1);  % Standard deviation for Vdmd

% % Case 4: companion DMD VCdmd
% mean_VCdmd = mean(abs(VCdmd), 1);  % Mean for VCdmd over 3880 trajectories
% std_VCdmd = std(abs(VCdmd), 0, 1);  % Standard deviation for VCdmd

% Case 5: sparsity-promoting DMD Humispdmd
mean_Vspdmd  = mean(real(Vspdmd ), 1);   % Mean for Vspdmd over 3880 trajectories
std_Vspdmd  = std(real(Vspdmd), 0, 1);  % Standard deviation for Vspdmd

% Plot Case 1: V0
% t= linspace(0,50, 51);
% Plot Case 1: Original data (mean_V0 and std_V0)
figure;
shadedErrorBar(t, mean_V0, std_V0, 'lineprops', {'--b', 'LineWidth', 1.5}, 'patchSaturation', 0.33);
hold on;

% % Plot Case 2: POD (mean_VPODTMs and std_VPODTMs) with orange color
% shadedErrorBar(t, mean_HPODTMs, std_HPODTMs, 'lineprops', {'--', 'Color', [1, 0.5, 0], 'LineWidth', 1.5}, 'patchSaturation', 0.075);  % Orange color [1, 0.5, 0]

% Plot Case 3: Vdmd (mean_Vdmd and std_Vdmd)
shadedErrorBar(t, mean_Vdmd, std_Vdmd, 'lineprops', {'-go', 'LineWidth', 1.5, 'MarkerFaceColor', 'g'}, 'patchSaturation', 0.33);

% Plot Case 4: VCdmd (mean_VCdmd and std_VCdmd)
% shadedErrorBar(t, mean_VCdmd, std_VCdmd, 'lineprops', {'-r', 'LineWidth', 1.5}, 'patchSaturation', 0.075);

% Plot Case 5: Vspdmd (mean_Humispdmd and std_Humispdmd)
shadedErrorBar(t, mean_Vspdmd, std_Vspdmd, 'lineprops', {'-mo', 'LineWidth', 1.5, 'MarkerFaceColor', 'w'}, 'patchSaturation', 0.33);

% Customize the plot with titles, labels, and a legend
title('Comparison of Temporal Dynamics');
xlabel('Time');
ylabel('Mean and Standard Deviation');

grid on;
legend({'Original', 'DMD', 'SPDMD'}, 'Location', 'Best');
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


%%
% Case 1: V0
mean_V0 = mean(real(V0), 1);   % Mean over all 3880 trajectories, for each time step
std_V0 = std(real(V0), 0, 1);  % Standard deviation over all 3880 trajectories, for each time step

% Case 2: Vdmd
mean_Vdmd = mean(real(Vdmd), 1);   % Mean for Vdmd over 3880 trajectories
std_Vdmd = std(real(Vdmd), 0, 1);  % Standard deviation for Vdmd

% Case 3: VCdmd
%mean_VCdmd = mean(real(VCdmd), 1);  % Mean for VCdmd over 3880 trajectories
%std_VCdmd = std(real(VCdmd), 0, 1);  % Standard deviation for VCdmd

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
%shadedErrorBar(t, mean_VCdmd, std_VCdmd, 'lineprops', '-r', 'patchSaturation', 0.33);

% Plot Case 4: Vspdmd
shadedErrorBar(t, mean_Vspdmd, std_Vspdmd, 'lineprops', {'-mo', 'MarkerFaceColor', 'm'}, 'patchSaturation', 0.33);

% Customizations
xlabel('Time');
ylabel('Amplitude');
title('Comparison of Four Trajectories (V0, Vdmd, VCdmd, Vspdmd) with Shaded Error Bars');
grid on;
legend({'V0', 'Vdmd', 'VCdmd', 'Vspdmd'}, 'Location', 'Best');
hold off;



% Plot Loss vs. Sparsity-promoting weights gamma (left y-axis)
figure;
yyaxis left;
semilogx(answer.gamma, answer.Ploss, 'bo-', 'LineWidth', 1, 'MarkerSize', 7);
%plot(answer.gamma, answer.Ploss, 'bo-', 'LineWidth', 1, 'MarkerSize', 7);
xlabel('\gamma', 'Interpreter', 'tex');
ylabel('Performance Loss (%)', 'Interpreter', 'tex');
set(gca, 'YColor', 'b'); % Set y-axis color for Loss plot

% Plot Nz vs. Sparsity-promoting weights gamma (right y-axis)
yyaxis right;
semilogx(answer.gamma, answer.Nz, 'rd--', 'LineWidth', 1, 'MarkerSize', 7);
%plot(answer.gamma, answer.Nz, 'rd--', 'LineWidth', 1, 'MarkerSize', 7);
ylabel('Number of Non-Zero Modes', 'Interpreter', 'tex');
set(gca, 'YColor', 'r'); % Set y-axis color for N_z plot

% Set font size and properties for labels
%xlab = xlabel('\gamma', 'Interpreter', 'tex');
xlabel('\gamma', 'Interpreter', 'tex');
%set(xlab, 'FontName', 'cmr10', 'FontSize', 26);
%set(gca, 'FontName', 'cmr10', 'FontSize', 20);

% Add a legend to distinguish between the two plots
legend('Loss', '# Modes', 'Location', 'best');
title('Accuracy vs. Model Reduction')
grid on;



%---------------
% Plot Loss vs. Sparsity-promoting weights gamma (left y-axis)
figure;
yyaxis left;
semilogx(answer.gamma, answer.Ploss, 'bo-', 'LineWidth', 1, 'MarkerSize', 7);
%plot(answer.gamma, answer.Ploss, 'bo-', 'LineWidth', 1, 'MarkerSize', 7);
xlabel('\gamma', 'Interpreter', 'tex');
ylabel('Performance Loss (%)', 'Interpreter', 'tex');
set(gca, 'YColor', 'b'); % Set y-axis color for Loss plot

% Plot Nz vs. Sparsity-promoting weights gamma (right y-axis)
yyaxis right;
semilogx(answer.gamma, answer.Nz, 'rd--', 'LineWidth', 1, 'MarkerSize', 7);
%plot(answer.gamma, answer.Nz, 'rd--', 'LineWidth', 1, 'MarkerSize', 7);
ylabel('Number of Non-Zero Modes', 'Interpreter', 'tex');
set(gca, 'YColor', 'r'); % Set y-axis color for N_z plot

% Set font size and properties for labels
%xlab = xlabel('\gamma', 'Interpreter', 'tex');
xlabel('\gamma', 'Interpreter', 'tex');
%set(xlab, 'FontName', 'cmr10', 'FontSize', 26);
%set(gca, 'FontName', 'cmr10', 'FontSize', 20);

% Add a legend to distinguish between the two plots
legend('Loss', '# Modes', 'Location', 'best');
title('Accuracy vs. Model Reduction')
grid on;

