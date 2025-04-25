%% superposotion of selected modes
% plot_superposition.m
% created by Z. ZHnag and Y. Susuki,

%% sort amplitudes for xpol, xsp, xdmd

[Norm_xdmd,Index_xdmd] = sort(abs(xdmd),'descend');
DEv_xdmd = Edmd(Index_xdmd);   %discrete-eigenvalues 
xdmd_sorted = xdmd(Index_xdmd);
DMDModes_xdmd = Phi(:,Index_xdmd);


kk =354; 
rr=answer.Nz(kk);

sort_xsp = answer.xsp(:,kk);  
%sort of the large amplitudes rather than rank of Eigvals
[Norm_xsp,Index_xsp] = sort(abs(answer.xsp(:,kk)),'descend');  %return the value of order/peak of amplitudes
DEv_xsp = Edmd(Index_xsp);   %discrete-eigenvalues 
xsp_sorted = answer.xsp(Index_xsp,kk);
DMDModes_xsp = Phi(:,Index_xsp);


[Norm_xpol,Index_xpol] = sort(abs(answer.xpol(:,kk)),'descend');
DEv_xpol = Edmd(Index_xpol);   %discrete-eigenvalues 
xpol_sorted = answer.xpol(Index_xpol,kk);
DMDModes_xpol = Phi(:,Index_xpol);


% give a specific spatial loctaion and related highlight index
% Define the row and column of the point you want to highlight, in[40,97]
% highlight_row = 4;  % <40 
 %highlight_col = 88;   % 97

 highlight_col = 4;  % <40 %ny
 highlight_row = 88;   % 97 % nz

%% option 1:
% Compute the linear indices for use in vectorized data
%highlight_index = (highlight_row - 1) * 97 + highlight_col;
% highlight_index_a = (highlight_row_a - 1) * 97 + highlight_col_a;
%highlight_index   = (highlight_col - 1) * 40 + highlight_row;
% highlight_index_a = (highlight_col_a - 1) * 40 + highlight_row_a;

%% Option2 === Highlight locations (in grid indices) ===

% highlight_pos     = [highlight_col, highlight_row];    % [z, y]
% highlight_pos_a   = [highlight_col_a, highlight_row_a];% [z, y]
% 
% % === Spatial grid settings ===
% ny = 40;  % y-direction (columns)
% nz = 97;  % z-direction (rows)
% 
% % Compute linear indices using sub2ind (used when flattening with reshape(..., [ny, nz])')
% highlight_index   = sub2ind([nz, ny], highlight_pos(1), highlight_pos(2));
% highlight_index_a = sub2ind([nz, ny], highlight_pos_a(1), highlight_pos_a(2));


%% Reconstraucted tempraol dynamics
 t = linspace(0, size(V0, 2), size(V0, 2));
time_dynamics = zeros(r, length(t));       % SPDMD reconstruted 
time_dynamics_Org = zeros(r, length(t));   % DMD reconstrauted 

% Compute the time dynamics for each time step
for iter = 1:length(t)
    time_dynamics(:, iter) = xsp_sorted.* exp(log(DEv_xsp) * t(iter));
    time_dynamics_Org(:, iter) = xdmd_sorted  .* exp(log(DEv_xdmd) * t(iter));
end

%% ------------ Compute Vspdmd---------------------------------
%Vspdmd = DMDModes_xsp * time_dynamics;  % sorted SPDMD 
Vspdmd = Phi*diag(answer.xsp(:,kk))*Vand;
Vspdmd_r_7 = Phi*diag(answer.xsp(:,400))*Vand;
Vspdmd_r_55 = Phi*diag(answer.xsp(:,324))*Vand;
Vspdmd_r_70 = Phi*diag(answer.xsp(:,231))*Vand;
Vspdmd_r_84 = Phi*diag(answer.xsp(:,1))*Vand;
%save('Vspdmd.mat','Vspdmd');

%% ----------- Compute Vdmd ------------------------------------
Vdmd = DMDModes_xdmd * time_dynamics_Org; %sorted DMD reconstructed
%Vdmd=Phi*diag(xdmd)*Vand;
%save('Vdmd.mat','Vdmd');

%% plot 
% % Calculate the linear index corresponding to this location
% highlight_index = (highlight_col - 1) * 40 + highlight_row;

% Create a single figure for combined plots
figure;
% Plot the temporal dynamics for the specific highlight index

plot(t, abs(V0(highlight_index, :)), 'LineWidth', 2, 'DisplayName', ['Original ${\bf{y}}_{k,{\rm{i}}}$, Index ', num2str(highlight_index)]);
hold on;
plot(t, abs(Vdmd(highlight_index, :)), 'LineWidth', 2, 'DisplayName', ['DMD ${\bf{y}}_{k,{\rm{i}}}$ ', num2str(highlight_index)]);
hold on;
plot(t, abs(Vspdmd_r_7(highlight_index, :)), 'LineWidth', 2, 'DisplayName', ['SPDMD-Superpoition y_k,i ', num2str(highlight_index)]);
hold on;
plot(t, abs(Vspdmd(highlight_index, :)), 'LineWidth', 2, 'DisplayName', ['SPDMD ${\bf{y}}_{k,{\rm{i}}}$ ', num2str(highlight_index)]);
hold on;
plot(t, abs(Vspdmd_r_55(highlight_index, :)), 'LineWidth', 2, 'DisplayName', ['SPDMD ${\bf{y}}_{k,{\rm{i}}}$', num2str(highlight_index)]);
hold on;
plot(t, abs(Vspdmd_r_70(highlight_index, :)), 'LineWidth', 2, 'DisplayName', ['SPDMD ${\bf{y}}_{k,{\rm{i}}}$', num2str(highlight_index)]);
hold on;
plot(t, abs(Vspdmd_r_84(highlight_index, :)), 'LineWidth', 2, 'DisplayName', ['SPDMD ${\bf{y}}_{k,{\rm{i}}}$ ', num2str(highlight_index)]);


% Customize the plot
title('Superposition Comparison');
xlabel('Time $k$');
ylabel('Evolution of $\sum_{j=1}^{r}{\rm{Re}}({[\phi}_j]^{\rm{i}}\lambda_j^kb_j)$');
legend('show'); % Automatically uses the 'DisplayName' property to create legend entries
grid on;
hold off;


%
%% optional: realx the optional
% Case 1: oiginal data V0
mean_V0 = mean(real(V0), 1);   % Mean over all 3880 trajectories, for each time step
std_V0 = std(real(V0), 0, 1);  % Standard deviation over all 3880 trajectories, for each time step

% Case 2: DMD reconstruction
mean_Vdmd = mean(real(Vdmd), 1);   % Mean for Vdmd over 3880 trajectories
std_Vdmd = std(real(Vdmd), 0, 1);  % Standard deviation for Vdmd

% % Case 3: companion DMD VCdmd
% mean_VCdmd = mean(real(VCdmd), 1);  % Mean for VCdmd over 3880 trajectories
% std_VCdmd = std(real(VCdmd), 0, 1);  % Standard deviation for VCdmd

% Case 4: sparsity-promoting DMD Humispdmd
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
