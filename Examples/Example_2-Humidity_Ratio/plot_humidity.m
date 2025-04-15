% plot for humidity data

% sort amplitudes for xpol, xsp, xdmd

[Norm_xdmd,Index_xdmd] = sort(abs(xdmd),'descend');
DEv_xdmd = Edmd(Index);   %discrete-eigenvalues 
xdmd_sorted = xdmd(Index_xdmd);  
DMDModes_xdmd = Phi(:,Index_xdmd);

kk=170;
rr=answer.Nz(kk);

%sort_xsp   
% sort of the large amplitudes rather than rank of Eigvals
[Norm_xsp,Index_xsp] = sort(abs(answer.xsp(:,kk)),'descend');  %return the value of order/peak ofamplitudes
DEv_xsp = Edmd(Index_xsp);   %discrete-eigenvalues 
xsp_sorted = answer.xsp(Index_xsp,kk);  % sorted the amplitudes
DMDModes_xsp = Phi(:,Index_xsp);


[Norm_xpol,Index_xpol] = sort(abs(answer.xpol(:,kk)),'descend');
DEv_xpol = Edmd(Index_xpol);   %discrete-eigenvalues 
xpol_sorted = answer.xpol(Index_xpol,kk);
DMDModes_xpol = Phi(:,Index_xpol);






% make a moive or demo for recosntructed spDMD
gif_filename = 'KM_spdmd_humidity.gif'; % Output GIF file name

% Define data dimensions
% rows = 40;
% cols = 97;

 rows =97;
 cols=40;

% Create and capture each frame
for i = 1:2:11
    mode_i = abs(DMDModes_xsp(:, i));
    mode_i = mode_i(1:cols*rows); % Adjust size if needed
    %reshap_mode_i= reshape(mode_i, [cols, rows]);
     reshap_mode_i= reshape(mode_i, [97, 40]);
 

    % Plot the spatial structure
     imagesc(data.y, data.z, reshap_mode_i); % humity
    colorbar; % Add colorbar for reference
   % axis xy; % Ensure the correct orientation
    xlim([0 2e4]);
    ylim([0 2e4]);
     xlabel("y");
    ylabel("z"); 
    title(['Humidity Spatial Mode ', num2str(i)]);
   % colormap('jet');
     colormap(brighten(redblueTecplot(21),-0.15));
    colorbar;
    % Adjust axis to fit data
   % axis equal; % Ensure aspect ratio is maintained

    
    % Capture the frame
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    % Write the frame to the GIF file
    if i == 1
        imwrite(imind, cm, gif_filename, 'gif', 'LoopCount', Inf, 'DelayTime', 1);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1);
    end
    
    % Close the figure
   % close(gcf);
end


%%% tempporal modes 
% Define modes to plot (example: 1, 3, 5, 7, 9, 11)
% kk =110;
% rr=answer.Nz(kk);
modesToPlot = 1:2:rr; % Adjust this based on your actual mode range
numModes = length(modesToPlot);

% Number of rows and columns in the subplot grid
numRows = 2;
numCols = 3;

% Define the row and column of the point you want to highlight
highlight_row = 31; % < 40 (this should be between 1 and 40)
highlight_col =90;  % < 97 (this should be between 1 and 97)

% Ensure the row and column are valid
% if highlight_row < 1 || highlight_row > 40 || highlight_col < 1 || highlight_col > 97
%     error('highlight_row must be between 1 and 40, and highlight_col must be between 1 and 97.');
% end

% Calculate the linear index corresponding to this location
highlight_index = (highlight_row - 1) * 97 + highlight_col; % 3880 size grid: 40 rows and 97 columns

% Create a figure for the subplots
figure;

for k = 1:numModes
    i = modesToPlot(k);
    
    % Determine subplot position
    subplotRow = ceil(k / numCols);
    subplotCol = mod(k - 1, numCols) + 1;
    
    % Extract and reshape mode_i
    % Here mode_i should be defined based on your data
    mode_i = abs(DMDModes_xdmd(:, i)); % Example: Replace with the correct mode extraction
    mode_i = mode_i(1:40 * 97); % Extract the first 40*97 elements
    
    % Plot the spatial structure in the subplot
    subplot(numRows, numCols, (subplotRow - 1) * numCols + subplotCol);
    load flujet % Load the data for the plotting
    imagesc(data.y, data.z, reshape(mode_i, [97, 40])); % Reshape and plot the spatial mode
    colorbar; % Add colorbar for reference
    colormap(brighten(redblueTecplot(110), -0.35)); % Adjust colormap
    xlim([0 2e4]);
    ylim([0 2e4]);
    xlabel("y");
    ylabel("z"); 
    title(['Spatial Mode ', num2str(i)]);
    
    % Highlight the specific point on the spatial mode
    hold on; % Keep the plot active
    highlight_y = data.y(highlight_row); % Get y-coordinate
    highlight_z = data.z(highlight_col); % Get z-coordinate
    plot(highlight_y, highlight_z, 'ko', 'MarkerSize', 9, 'MarkerFaceColor', 'yellow'); % Use a yellow circle for highlighting
    box;
    hold off; % Release the plot
end



%% Plot temporal modes 

t= linspace(0,49, 50);
time_dynamics = zeros(r,length(t));
%sort_time_dynamics = zeros(r,length(t));
for iter = 1:length(t)
%time_dynamics (:,iter) = (xdmd_sorted.*exp(log(DEv_xdmd)*t(iter))); original
%time_dynamics (:,iter) = (xpol_sorted.*exp(log(DEv_xpol)*t(iter)));
%time_dynamics (:,iter) = (xsp_sorted.*exp(log(DEv_xsp)*t(iter)));
% time_dynamics_Org(:, iter) = answer.xsp(:,kk).* exp(log(DEv_xdmd) * t(iter)); % without sorted 
 time_dynamics(:, iter) =  xsp_sorted.* exp(log(DEv_xsp) * t(iter)); % sorted
end

% Humispdmd = DMDModes_xsp * time_dynamics;  % sorted 
  Humispdmd = Phi * diag(answer.xsp(:,kk))*Vand; % withoutsoretd
% save data 
save('Humispdmd',"Humispdmd");


% Visualize Temporal Dynamics 
% for i = 1:2:rr
%     figure; % Create a new figure for each mode
%     plot(t, real(time_dynamics(i, :))); % Plot real part of temporal dynamics
%     title(['Temporal Dynamics of Mode ', num2str(i)]);
%     xlabel('Time');
%     ylabel('Evolution $\lambda_j^{k}b_{j}$');
% end

%% Option 1: plot the temperal on one figure
% Define line styles and markers for each mode

line_styles = {'-', '--', ':', '-.', '-', '--'};
markers = {'o', 's', 'd', '^', 'v', 'x'};

% Plot Temporal Dynamics for each mode on the same figure
figure;
hold on; % Keep all plots in the same figure

for i = 1:2:rr
    plot(t, real(time_dynamics(i, :)), 'LineStyle', line_styles{mod(i, length(line_styles))+1}, ...
        'Marker', markers{mod(i, length(markers))+1}, 'LineWidth', 1.5, 'DisplayName', ['Mode ', num2str(i)]);
end

title('Temporal Dynamics of Selected Modes');
xlabel('Time');
ylabel('Evolution $\lambda_j^{k}b_{j}$', 'Interpreter', 'latex');
legend('show'); % Display the legend
box;
hold off; % Release the figure

title('Temporal Dynamics of Selected Modes');
xlabel('Time');
ylabel('Evolution $\lambda_j^{k}b_{j}$', 'Interpreter', 'latex');
legend('show'); % Display the legend
box on; % Show the box around the plot
hold off; % Release the figure


%% Option 2: Compared with diverse methods 

%============DMD=======================
% t= linspace(0,50, 51);
time_dynamics_Org = zeros(r, length(t));

% Compute the time dynamics for each time step
for iter = 1:length(t)
  %     time_dynamics_Org(:, iter) = xdmd  .* exp(log(Edmd) * t(iter));
 %time_dynamics_Org(:, iter) = Norm_xdmd.*exp(log(DEv_xdmd)*t(iter)); original
 time_dynamics_Org(:, iter) = xdmd.*exp(log(DEv_xdmd)*t(iter)); original
end

% Compute Humidity ratio via DMD 
Humidmd = Phi * diag(xdmd)*Vand;;
%Humidmd =  DMDModes_xdmd * time_dynamics_Org;

% save data
save('Humidmd',"Humidmd");

%======Companion DMD============================
% Truncated data (snapshots) from Hdata
 Humfull = Hdata(:,1:51);
 [HCKModes,HCKEv,HCNorms,HCdmd] = CompanionMatrix_DMD( Humfull);


%======Classic POD============================
[HPODMs, HPODEvals, HSortPOD, HPODTMs] = classicPOD(Humfull);



%===== Residual DMD==========================
% %[G,K,L,PX,PY,PSI_x,PSI_y,PSI_y2] = kernel_ResDMD(V0,V1,'type',"Gaussian");
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



%% Plot temporal dynamics (specific laocation and mean)
% Plot the temporal dynamics for the specific highlight index
figure;
plot(t, real(V0(highlight_index, :)), 'LineWidth', 2, 'DisplayName', ['Original Index ', num2str(highlight_index)]);
hold on;

plot(t, real(HPODTMs(highlight_index, 1:end-1)), 'LineWidth', 2, 'DisplayName', ['SPDMD-Recon y_k,i ', num2str(highlight_index)]);
hold on;

plot(t, real(HCdmd(highlight_index, :)), 'LineWidth', 2, 'DisplayName', ['CDMD-Recon y_k,i ', num2str(highlight_index)]);
hold on;

plot(t, real(Humidmd(highlight_index, :)), 'LineWidth', 2, 'DisplayName', ['DMD-Recon y_k,i ', num2str(highlight_index)]);
hold on;

plot(t, real(Humispdmd(highlight_index, :)), 'LineWidth', 2, 'DisplayName', ['SPDMD-Recon y_k,i ', num2str(highlight_index)]);
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
% rmse_Humidmd = sqrt(mean((real(V0(highlight_index, :) - Humidmd(highlight_index, :))).^2, 1))/ sqrt(length(t));
% 
% % Compute RMSE over time for SPDMD-reconstructed data
% rmse_Humispdmd = sqrt(mean((real(V0(highlight_index, :) - Humispdmd(highlight_index, :))).^2, 1))/ sqrt(length(t));
% 
% % Plot the RMSE over time for DMD and SPDMD
% figure;
% plot(t, rmse_Humidmd, 'LineWidth', 2, 'DisplayName', 'RMSE DMD');
% hold on;
% plot(t, rmse_Humispdmd, 'LineWidth', 2, 'DisplayName', 'RMSE SPDMD');
% xlabel('Time');
% ylabel('RMSE');
% title(['Root Mean Square Error ', num2str(highlight_index)]);
% legend('show');
% grid on;
% hold off;



% Case 1: oiginal Humidity ratio data V0
mean_V0 = mean(abs(V0), 1);   % Mean over all 3880 trajectories, for each time step
std_V0 = std(abs(V0), 0, 1);  % Standard deviation over all 3880 trajectories, for each time step

% Case 2: POD reconstruction
mean_HPODTMs = mean(abs(HPODTMs(:,1:end-1)), 1);   % Mean for Vdmd over 3880 trajectories
std_HPODTMs = std(abs(HPODTMs(:,1:end-1)), 0, 1);  % Standard deviation for Vdmd

% Case 3: DMD reconstruction
mean_Humidmd = mean(abs(Humidmd), 1);   % Mean for Vdmd over 3880 trajectories
std_Humidmd = std(abs(Humidmd), 0, 1);  % Standard deviation for Vdmd

% Case 4: companion DMD HCdmd
mean_HCdmd = mean(abs(HCdmd), 1);  % Mean for VCdmd over 3880 trajectories
std_HCdmd = std(abs(HCdmd), 0, 1);  % Standard deviation for VCdmd

% Case 5: sparsity-promoting DMD Humispdmd
mean_Humispdmd  = mean(abs(Humispdmd ), 1);   % Mean for Vspdmd over 3880 trajectories
std_Humispdmd  = std(abs(Humispdmd), 0, 1);  % Standard deviation for Vspdmd

% Plot Case 1: V0
% t= linspace(0,50, 51);
% Plot Case 1: Original data (mean_V0 and std_V0)
figure;
shadedErrorBar(t, mean_V0, std_V0, 'lineprops', {'--b', 'LineWidth', 1.5}, 'patchSaturation', 0.33);
hold on;

% Plot Case 2: POD (mean_HPODTMs and std_HPODTMs) with orange color
shadedErrorBar(t, mean_HPODTMs, std_HPODTMs, 'lineprops', {'--', 'Color', [1, 0.5, 0], 'LineWidth', 1.5}, 'patchSaturation', 0.075);  % Orange color [1, 0.5, 0]

% Plot Case 3: Hdmd (mean_Humidmd and std_Humidmd)
shadedErrorBar(t, mean_Humidmd, std_Humidmd, 'lineprops', {'-go', 'LineWidth', 1.5, 'MarkerFaceColor', 'g'}, 'patchSaturation', 0.33);

% Plot Case 4: HCdmd (mean_HCdmd and std_HCdmd)
shadedErrorBar(t, mean_HCdmd, std_HCdmd, 'lineprops', {'-r', 'LineWidth', 1.5}, 'patchSaturation', 0.075);

% Plot Case 5: Hspdmd (mean_Humispdmd and std_Humispdmd)
shadedErrorBar(t, mean_Humispdmd, std_Humispdmd, 'lineprops', {'-mo', 'LineWidth', 1.5, 'MarkerFaceColor', 'w'}, 'patchSaturation', 0.33);

% Customize the plot with titles, labels, and a legend
title('Comparison of Temporal Dynamics');
xlabel('Time');
ylabel('Mean and Standard Deviation');

grid on;
legend({'Original', 'POD', 'DMD', 'Companion DMD', 'SPDMD'}, 'Location', 'Best');
hold off;


%% sort or account the number of zero amplitude related eigenvals 
 nonzero_indices = find(answer.xsp(:,kk)); % Step 1: Find nonzero amplitudes
  [~, sort_order] = sort(nonzero_indices, 'descend'); 
  ival = nonzero_indices(sort_order);    % also ralted to `DEv_xsp' (the order of eigvals)

 
%% plot eigenvalues vs spdmd
figure;
% % Plot POD eigenvalues (HPODEvals) with black diamonds
% plot(real(HPODEvals), imag(HPODEvals), 'ks', 'MarkerSize', 10, 'LineWidth', 2);
% hold on;

% Plot DMD eigenvalues (Edmd) with blue circles
plot(real(Edmd), imag(Edmd), 'bo', 'MarkerSize', 8, 'LineWidth', 1.5);
hold on;

% Plot HCK eigenvalues (HCKEv) with black squares
plot(real(HCKEv), imag(HCKEv), 'kd', 'MarkerSize', 8, 'LineWidth', 1.5);
hold on;

% Highlight selected DMD eigenvalues (Edmd(ival)) with red crosses
plot(real(Edmd(ival)), imag(Edmd(ival)), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
hold on;

% Draw a unit circle
rectangle('Position', [-1 -1 2 2], 'Curvature', 1, 'EdgeColor', 'k', 'LineStyle', '--');

% Customize axis labels and title
xlabel('Real part');
ylabel('Imaginary part');
title('Eigenvalues Humidity Ratio');

% Set axis limits and aspect ratio
axis(1.2*[-1 1 -1 1]);
axis square;
grid on;




%% plot growth rate 
% Compute growth rates from eigenvalues
growth_rates = real(log(Edmd)); % Growth rate (real part of log eigenvalue)

% Define sparsity-level growth rates
sp_growth_rates = real(log(Edmd(ival))); % Specific modes

% Plot growth rates
figure;
hold on;

% Plot all modes and sparsity-level modes
plot(growth_rates, imag(log(Edmd)), 'bo', ... % All modes
     sp_growth_rates, imag(log(Edmd(ival))), 'r+', ... % Sparsity-level modes
     'LineWidth', 1, 'MarkerSize', 7);

% Add vertical line at x = 0
line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--');

% Add horizontal line at y = 0
line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--');

% Shade the region where x <= 0
x = [min(xlim), 0, 0, min(xlim)];
y = [min(ylim), min(ylim), max(ylim), max(ylim)];
fill(x, y, [0.9, 0.9, 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Labels and formatting
xlabel('Growth Rate', 'interpreter', 'tex');
ylabel('Frequency (Imaginary Part of log-Eigenvalues)', 'interpreter', 'tex');
set(gca, 'FontName', 'cmr10', 'FontSize', 20);
legend('DMD', 'spDMD');
hold off;



%% --show Norm-e-folding time periodic

% Assuming eigenvalues are given in Edmd(ival)
eigenvalues = Edmd(ival); % Replace 'ival' with the index range of interest
                           % mainly from DEv_xsp 

% Preallocate arrays for results
num_eigenvalues = length(eigenvalues);
norms = zeros(num_eigenvalues, 1);
e_folding_times = zeros(num_eigenvalues, 1);
periods = zeros(num_eigenvalues, 1);

% Define constants
time_step = 1; % Adjust this based on the time step of your system (e.g., days, months)

% Loop through each eigenvalue
for i = 1:num_eigenvalues
    lambda = DEv_xsp(i);
    
    % Compute norm
    norms(i) = abs(lambda);
    
    % Compute e-folding time
    if abs(lambda) < 1
        e_folding_times(i) = -time_step / log(abs(lambda)); % Only compute for stable modes
    else
        e_folding_times(i) = Inf; % Unstable modes do not have an e-folding time
    end
    
    % Compute period (if applicable)
    omega = imag(log(lambda)); % Imaginary part of the logarithm gives angular frequency
    if omega ~= 0
        periods(i) = 2 * pi / abs(omega); % Period from angular frequency
    else
        periods(i) = Inf; % Non-oscillatory modes do not have a period
    end
end


%% show the print results 
% Display results with appropriate formatting
fprintf(' Mode  Index_Mode  Norm_xsp  Norm_Eigvals   DEv_xsp   E-Folding Time       Period\n');
fprintf('------------------------------------------------------------------------------------------\n');

for i = 1:num_eigenvalues
    % Split into separate lines for better clarity
    fprintf('%5d      %5d      %.2f       %.2f      %.3f + %.3fi', ...
            i, Index_xsp(i), Norm_xsp(i), abs(DEv_xsp(i)), real(DEv_xsp(i)), imag(DEv_xsp(i)));
    
    fprintf('       %.2f                %.2f\n', ...
            e_folding_times(i), periods(i));
end



%%  Make Cross verification

figure;

% Plot Loss vs. Sparsity-promoting weights gamma (left y-axis)
yyaxis left;
%semilogx(answer.gamma, answer.Ploss, 'bo-', 'LineWidth', 1, 'MarkerSize', 7);
plot(answer.gamma, answer.Ploss, 'bo-', 'LineWidth', 1, 'MarkerSize', 7);
xlabel('Sparsity-promoting weights \gamma', 'Interpreter', 'tex');
ylabel('Performance Loss (%)', 'Interpreter', 'tex');
set(gca, 'YColor', 'b'); % Set y-axis color for Loss plot

% Plot Nz vs. Sparsity-promoting weights gamma (right y-axis)
yyaxis right;
semilogx(answer.gamma, answer.Nz, 'rd--', 'LineWidth', 1, 'MarkerSize', 7);
ylabel('N_z (Number of Non-Zero Modes)', 'Interpreter', 'tex');
set(gca, 'YColor', 'r'); % Set y-axis color for N_z plot

% Set font size and properties for labels
%xlab = xlabel('\gamma', 'Interpreter', 'tex');
xlabel('\gamma', 'Interpreter', 'latex');
set(gca, 'FontName', 'cmr10', 'FontSize', 20);

% Add a legend to distinguish between the two plots
legend('Loss', '# Modes', 'Location', 'best');
title('Accuarcy vs. Model Reduction')
grid on;

