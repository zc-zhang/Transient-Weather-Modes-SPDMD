% plot_spdmd_scalar.m


kk =50; %vorticity
%kk=123;
rr=answer.Nz(kk);
% figure
% plot(abs(xdmd),'LineWidth',2.5,'Color','b','LineStyle','-')
% hold on 
% %kk=243;
% plot(abs(answer.xpol(:,kk)),'LineWidth',2.5,'Color','r','LineStyle',':')
% %plot(abs(answer.xsp(:,224)),'LineWidth',2.5,'Color','g','LineStyle',':')
% %legend('DMD mode amplitude');  
% ylabel('Amplitudes of abs(b) vs abs(bsp)');
% xlabel('Mode number');

% Spectrum of DT/CT for DMD
% figure('Position', [100 100 600 300]);
% subplot(1,2,1);
% plot(Edmd , 'rx');
% rectangle('Position', [-1 -1 2 2], 'Curvature', 1, ...
% 'EdgeColor', 'b', 'LineStyle', '--');
% xlabel('Real part');
%  ylabel('Imaginary part');
%  title('Discrete Eigenvalues of DMD Edmd')
% axis (1.2*[-1 1 -1 1]);
% axis square;
% grid on
% % CT spectrum
% subplot(1,2,2);
% plot(omega, 'k.');
% line ([0 0], 200*[-1 1], 'Color', 'b', 'LineStyle', '--');
%  xlabel('Real part');
%  ylabel('Imaginary part');
% title('Continuous Eigenvalues of DMD omega');
% axis ([-5 1 -5 +5]);
% axis square;
% grid on


%figure
% xdmd_sp = answer.xsp;
% xdmd_pol =answer.xpol;
% plot(abs(xdmd),'LineWidth',2,'Color','k','LineStyle','-')
% %legend('DMD mode amplitude xdmd','--');  
% hold on
% plot(abs(xdmd_sp),'LineWidth',2,'Color','b','LineStyle','-.')
% %legend('DMD mode amplitude xdmdsp',':');  
% hold on
% plot(abs(xdmd_pol),'LineWidth',2,'Color','r','LineStyle',':')
% %legend('DMD mode amplitude xdmdpol');  
% ylabel('abs(xdmd)');
% xlabel('Mode number');

% |xdmd| vs frequency
figure;
plot(imag(log(Edmd)),abs(xdmd),'bo')         % log(Edmd)/dt = omega, dt=1
xlab = xlabel('frequency','interpreter','tex')
set(xlab,'FontName','cmmr10','FontSize',26)
ylab = ylabel('amplitude','interpreter','tex')
set(ylab,'FontName','cmmr10','FontSize',26)
h = get(gcf,'CurrentAxes'); 
set(h,'FontName','cmr10','FontSize',20,'xscale','lin','yscale','lin')

% |xdmd| vs real part
figure;
plot(real(log(Edmd)),abs(xdmd),'ro')
xlab = xlabel('real','interpreter','tex')
set(xlab,'FontName','cmmr10','FontSize',26)
ylab = ylabel('amplitude','interpreter','tex')
set(ylab,'FontName','cmmr10','FontSize',26)
h = get(gcf,'CurrentAxes'); 
set(h,'FontName','cmr10','FontSize',20,'xscale','lin','yscale','lin')

% Performance loss for the polished vector of amplitudes vs gamma 
figure;
semilogx(answer.gamma,answer.Ploss,'ko','LineWidth',1,'MarkerSize',7)
xlab = xlabel('\gamma','interpreter','tex')
set(xlab,'FontName','cmr10','FontSize',26)
ylab = ylabel('performance loss (%)','interpreter','tex')
set(ylab,'FontName','cmr10','FontSize',26)
h = get(gcf,'CurrentAxes'); 
set(h,'FontName','cmr10','FontSize',20)
axis([answer.gamma(1) answer.gamma(end) 0 1.05*answer.Ploss(end)])

% Number of non-zero amplitudes vs gamma
figure;
semilogx(answer.gamma,answer.Nz,'rd','LineWidth',1,'MarkerSize',7)
xlab = xlabel('\gamma','interpreter','tex')
set(xlab,'FontName','cmr10','FontSize',26)
ylab = ylabel('N_z','interpreter','tex')
set(ylab,'FontName','cmr10','FontSize',26)
h = get(gcf,'CurrentAxes'); 
set(h,'FontName','cmr10','FontSize',20)
%axis([answer.gamma(1) answer.gamma(end) 0 1.05*answer.Nz(1)])


% the no of non-zero amplitues-cardnality card(xsp) vs vairous gamma level
figure;
plot(answer.Nz,answer.Nz,'bo')
xlabel('Sparsity-promoting weights \gamma','interpreter','tex')
ylabel('Cardnality card(b)','interpreter','tex')

figure;
%plot(answer.gamma,answer.Ploss,'bo')
semilogx(answer.gamma,answer.Ploss,'bo')
xlabel('Sparsity-promoting weights \gamma','interpreter','tex')
ylabel('Loss','interpreter','tex')


% Spectrum of DT system for a certain value of gamma
%answer.Nz(kk) % index of non-zero amplitudes in gamma 
%xsp = answer(xsp); 
%==================================
nonzero_indices = find(answer.xsp(:,kk)); % Step 1: Find nonzero amplitudes
[~, sort_order] = sort(nonzero_indices, 'descend'); 
ival = nonzero_indices(sort_order);    % also ralted to `DEv_xsp' (the order of eigvals)
%=================================
%ival = find(answer.xsp(:,kk));
% ival = Index_xsp(1:rr); % according to Index_xsp to show the orders of ampli
% plot(real(Edmd),imag(Edmd),'ko',real(Edmd(ival)),imag(Edmd(ival)),'r+', ...
%     'LineWidth',1,'MarkerSize',15)
% xlab = xlabel('Re($\lambda_i$)','interpreter','latex')
% set(xlab,'FontName','cmr10','FontSize',26)
% ylab = ylabel('Im($\lambda_i$)','interpreter','latex')
% set(ylab,'FontName','cmr10','FontSize',26)
% %h = get(gcf,'CurrentAxes'); 
% set(h,'FontName','cmr10','FontSize',10,'xscale','lin','yscale','lin')
% hold
% % plot a unit circlek
% format compact                    % tighten loose format
% format long e                     % make numerical output in double precision
% theta = linspace(0,2*pi,100);     % create vector theta
% x = cos(theta);                   % generate x-coordinate
% y = sin(theta);                   % generate y-coordinate
% plot(x,y,'b--','LineWidth',2);    % plot unit circle
% axis('equal');
% hold

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

%% Eigenvalues using scatter 

answer.Nz(kk)  % index of non-zero amplitudes in gamma 

% Use Index_xsp to show the order of amplitudes
ival = Index_xsp(1:rr);  % Select the first 'rr' indices according to sorted amplitudes

% Create a figure for plotting
figure;

% Scatter plot for all DMD eigenvalues (Edmd)
scatter(real(Edmd), imag(Edmd), 65, abs(Edmd), 'filled');  % '50' sets the marker size, 'abs(Edmd)' provides color by magnitude
hold on;

% Scatter plot for selected SpDMD eigenvalues (marked with stars)
scatter(real(Edmd(ival)), imag(Edmd(ival)), 100, 'r', 'p', 'LineWidth', 2);  % 'p' for star marker, red color

% Plot unit circle
rectangle('Position', [-1 -1 2 2], 'Curvature', 1, 'EdgeColor', 'k', 'LineStyle', '--');  % Unit circle border

% Add labels and title with formatting
xlabel('Re($\lambda_i$)', 'Interpreter', 'latex', 'FontSize', 30);
ylabel('Im($\lambda_i$)', 'Interpreter', 'latex', 'FontSize', 30);
title(' Koopman Eigenvalues SST', 'FontSize', 16);

% Set axis limits and make the axes square
axis(1.5*[-1 1 -1 1]);
axis square;

% Add grid lines
grid on;
box;
% Add colorbar to show the magnitude scale
colorbar;
colormap Sky;  % Set the colormap for better visualization (can be customized)

%% Growth rate
% Compute growth rates from eigenvalues
growth_rates = real(log(Edmd)); % Growth rate (real part of log eigenvalue)

% Define sparsity-level growth rates
sp_growth_rates = real(log(Edmd(ival))); % Specific modes

% Plot growth rates
figure;
plot(imag(log(Edmd)), growth_rates, 'ko', ... % All modes
     imag(log(Edmd(ival))), sp_growth_rates, 'r+', ... % Sparsity-level modes
     'LineWidth', 1, 'MarkerSize', 7);
xlabel('Frequency (Imag Part of log-Eigvals)', 'interpreter', 'tex');
ylabel('Growth Rate ', 'interpreter', 'tex');
set(gca, 'FontName', 'cmr10', 'FontSize', 20);
legend('DMD', 'spDMD');

%%% change the x andY label
growth_rates = real(log(Edmd)); % Growth rate (real part of log eigenvalue)

% Define sparsity-level growth rates
sp_growth_rates = real(log(Edmd(ival))); % Specific modes

% Plot growth rates  % use imag(log(Edmd)) or abs(log(Edmd)) or angle(log(Edmd))
figure;
plot(growth_rates, imag(log(Edmd)),  'ko', ... % All modes
     sp_growth_rates, imag(log(Edmd(ival))), 'r+', ... % Sparsity-level modes
     'LineWidth', 1, 'MarkerSize', 7);
xlabel('Growth Rate ', 'interpreter', 'tex');
ylabel('Frequency (Imag Part of log-Eigvals)', 'interpreter', 'tex');
set(gca, 'FontName', 'cmr10', 'FontSize', 20);
legend('DMD', 'spDMD');


%%%%%% Modifty for growth rate 

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

% 

% |xdmd| and |xpol| vs frequency for a certain value of gamma
% amplitudes in log scale 
answer.Nz(kk) % number of non-zero amplitudes
%ival = find(answer.xsp(:,kk));
figure;
%semilogy(imag(log(Edmd)),abs(xdmd),'ko', ...
    %  imag(log(Edmd(ival))),abs(answer.xpol(ival,kk)),'r+', ...
    % 'LineWidth',1,'MarkerSize',7)
 plot(imag(log(Edmd)),abs(xdmd),'ko', ...
      imag(log(Edmd(ival))),abs(answer.xpol(ival,kk)),'r+', ...
     'LineWidth',1,'MarkerSize',7)
xlab = xlabel('frequency','interpreter','tex')
set(xlab,'FontName','cmr10','FontSize',26)
ylab = ylabel('amplitude','interpreter','tex')
set(ylab,'FontName','cmr10','FontSize',26)
h = get(gcf,'CurrentAxes'); 
set(h,'FontName','cmr10','FontSize',20)


%%%%%%%% Growth/decay rate
% Compute growth/decay rates
growth_decay_rate = abs(Edmd);
sp_growth_decay_rate = abs(Edmd(ival));
% Plot the growth/decay rates
figure;
plot(imag(log(Edmd)), growth_decay_rate, 'ko', ...
     imag(log(Edmd(ival))), sp_growth_decay_rate, 'r+', ...
     'LineWidth', 1, 'MarkerSize', 7);
xlabel('Frequency','interpreter','tex');
ylabel('Growth/Decay Rate','interpreter','tex');
set(gca, 'FontName', 'cmr10', 'FontSize', 20);





% Determine data for performance loss vs number of dmd modes plots
Nz(1) = answer.Nz(1);
Ploss(1) = answer.Ploss(1);
ind = 1;

for i = 1:length(answer.gamma)-1,

    if (answer.Nz(i) == answer.Nz(i+1)),

        ind = ind;

    else 

        ind = ind+1;
        Nz(ind) = answer.Nz(i+1);
        Ploss(ind) = answer.Ploss(i+1);

    end


end

clear ind

% Performance loss vs number of dmd modes
figure
Nz = answer.Nz;
Ploss =answer.Ploss;
plot(Nz,Ploss,'ko')
xlab = xlabel('number of dmd modes','interpreter','tex')
set(xlab,'FontName','cmr10','FontSize',26)
ylab = ylabel('performance loss (%)','interpreter','tex')
set(ylab,'FontName','cmr10','FontSize',26)
h = get(gcf,'CurrentAxes'); 
set(h,'FontName','cmr10','FontSize',20)
axis([Nz(end) Nz(1) 0 1.05*Ploss(end)])


%% =============Seshape MODES of Ydmd ===============
% % Visualize Spatial Modes

%rr=25;
% for i = 1:rr
%     figure; % Create a new figure for each mode
%     % mode_i = abs(DMDModes_xdmd(:, i)); 
%     % mode_i = abs(DMDModes_xpol(:, i)); % Take absolute value to visualize complex values
%      mode_i = abs(DMDModes_xsp(:, i)); 
%     mode_i = mode_i(1:40*97); % Extract the first 10*60 elements (assuming each tmpV 40*97)
%     imagesc(reshape(mode_i, [40, 97])); % Plot spatial structure
%     title(['Spatial Structure of Mode ', num2str(i)]);
%     colormap('jet');
%     colorbar;
% end

% spDMD modes w.r.t. the order of eigvalues `DEv_xsp'
for i = 1:2:rr
    figure; % Create a new figure for each mode
    % mode_i = abs(DMDModes_xdmd(:, i)); 
    % mode_i = abs(DMDModes_xpol(:, i)); % Take absolute value to visualize complex values
     mode_i = abs(DMDModes_xsp(:, i)); 
    mode_i = mode_i(1:40*97); % Extract the first 10*60 elements (assuming each tmpV 40*97)
 imagesc(data.y, data.z, reshape(mode_i, [40, 97])');  % very important use transpose: vortex 
   % colormap(jet); % Choose a colormap
    colorbar; % Add colorbar for reference
     colormap(brighten(redblueTecplot(21),-0.55));
    axis xy; % Ensure the correct orientation
    xlim([0 2e4]);
    ylim([0 2e4]);
     xlabel("y");
    ylabel("z"); 
    title(['Spatial Structure of Mode ', num2str(i)]);
   % colorbar;
end


%----------------------------Subplot the above Spatial mode==========
% Define modes to plot (example: 1, 3, 5, 7, 9, 11)
% modesToPlot = 1:2:rr-1 % Adjust this based on your actual mode range
% numModes = length(modesToPlot);
% 
% % Number of rows and columns in the subplot grid
% numRows = 2;
% numCols = 3;
% 
% % Create a new figure for all subplots
% figure;
% 
% % Loop through the modes you want to plot
% for idx = 1:min(numModes, numRows * numCols) % Ensure we do not exceed subplot grid size
%     % Get the mode index
%     modeIndex = modesToPlot(idx);
% 
%     % Determine subplot position
%     subplotRow = ceil(idx / numCols); % Calculate row number
%     subplotCol = mod(idx - 1, numCols) + 1; % Calculate column number
% 
%     % Select the subplot position
%     subplot(numRows, numCols, (subplotRow - 1) * numCols + subplotCol);
% 
%     % Extract and process mode data
%     mode_i = abs(DMDModes_xsp(:, modeIndex)); 
%     mode_i = mode_i(1:40*97); % Extract the first 40*97 elements (assuming each tmpV 40*97)
% 
%     % Plot the data
%     imagesc(data.y, data.z, reshape(mode_i, [40, 97])'); % Transpose data
%     colorbar; % Add colorbar for reference
%     colormap(brighten(redblueTecplot(21), -0.55));
%     axis xy; % Ensure the correct orientation
%     xlim([0 2e4]);
%     ylim([0 2e4]);
%     xlabel("y");
%     ylabel("z");
%     title(['Spatial Structure of Mode ', num2str(modeIndex)]);
% end
%------------------------------------------------------------------------


t= linspace(0,42, size(V0, 2));
time_dynamics = zeros(r,length(t));
%sort_time_dynamics = zeros(r,length(t));
for iter = 1:length(t)
%time_dynamics (:,iter) = (Norm_xdmd.*exp(log(DEv_xdmd)*t(iter))); original
%time_dynamics (:,iter) = (Norm_xpol.*exp(log(DEv_xpol)*t(iter)));
time_dynamics (:,iter) = (Norm_xsp.*exp(log(DEv_xsp)*t(iter)));
%sort_time_dynamics(:,iter) = (b_sort.*exp(omega*t(iter)));
end

% Visualize Temporal Dynamics 
for i = 1:2:rr
    figure; % Create a new figure for each mode
    plot(t, real(time_dynamics(i, :))); % Plot real part of temporal dynamics
    title(['Temporal Dynamics of Mode ', num2str(i)]);
    xlabel('Time');
    ylabel('Evolution $\lambda_j^{k}b_{j}$');
end


%%%%% plot the seleted temoral in one figure

t = linspace(0, 42, size(V0, 2));
time_dynamics = zeros(r, length(t));

% Compute the time dynamics for each time step
for iter = 1:length(t)
    time_dynamics(:, iter) = Norm_xsp .* exp(log(DEv_xsp) * t(iter));
end

%%% DMDModes_xsp is the DMD modes (i.e., Phi) with sparsity case
Vspdmd = DMDModes_xsp*time_dynamics;


%% hightindex 
plot(real(Vspdmd(highlight_index,:)))
hold on 
% Visualize Temporal Dynamics (e.g., real, image, phase, angle...)
figure; % Create a new figure for plotting
hold on; % Retain plots so that new plots are added to the existing figure

for i = 1:2:rr
    plot(t, real(time_dynamics(i, :))); % Plot real part (or imag,abs,angle) of temporal dynamics for each mode
end

title('Temporal Dynamics of Modes');
xlabel('Time');
ylabel('Evolution $\lambda_j^{k}b_{j}$');
legend(arrayfun(@(x) ['Temporal Mode ', num2str(x)], 1:2:rr, 'UniformOutput', false)); % Add a legend for clarity
hold off; % Release the hold on the figure


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%%% temporal modes of selsted + the i-the element of recosntraucted Vspdmd
% Define time vector and initialize time dynamics
t = linspace(0, 42, size(V0, 2));
time_dynamics = zeros(r, length(t));

% Compute the time dynamics for each time step
for iter = 1:length(t)
    time_dynamics(:, iter) = Norm_xsp .* exp(log(DEv_xsp) * t(iter));
end

% Compute Vspdmd
Vspdmd = DMDModes_xsp * time_dynamics;

% Create a single figure for combined plots
figure;

% Plot the temporal dynamics for the specific highlight index
plot(t, real(Vspdmd(highlight_index, :)), 'LineWidth', 2, 'DisplayName', ['Highlight Index ', num2str(highlight_index)]);
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



%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% for moive the Koopman modes sparsity


% Parameters
%rr = 25; % Number of modes
% gif_filename = 'KM_spdmd_scale_magnitude01.gif'; % Output GIF file name
% 
% % Define data dimensions
% rows = 40;
% cols = 97;
% 
% % Create and capture each frame
% for i = 1:2:rr
%     mode_i = abs(DMDModes_xsp(:, i));
%     mode_i = mode_i(1:rows*cols); % Adjust size if needed
%     %reshaped_mode = reshape(mode_i, [rows, cols]);
%     % reshaped_mode_v = reshape(mode_i, [rows, cols]);
%     % reshaped_mode_v = reshape(mode_i, [rows, cols]);
% 
% 
%     % Plot the spatial structure
%     imagesc(data.y, data.z, reshape(mode_i, [rows, cols])'); % vortex
%      %colormap(jet); % Choose a colormap
%     colorbar; % Add colorbar for reference
%     axis xy; % Ensure the correct orientation
%     xlim([0 2e4]);
%     ylim([0 2e4]);
%      xlabel("y");
%     ylabel("z"); 
%     title(['Spatial Mode ', num2str(i)]);
%     %colormap('jet');
%       colormap(brighten(redblueTecplot(21),-0.55));
%     colorbar;
% 
%     % Adjust axis to fit data
%    % axis equal; % Ensure aspect ratio is maintained
%     %axis([1 cols 1 rows]); % Ensure axes fit the data dimensions
%     %xlim([1 97]);
%     %ylim([1 40]);
% 
%     % Capture the frame
%     frame = getframe(gcf);
%     im = frame2im(frame);
%     [imind, cm] = rgb2ind(im, 256);
% 
%     % Write the frame to the GIF file
%     if i == 1
%         imwrite(imind, cm, gif_filename, 'gif', 'LoopCount', Inf, 'DelayTime', 1);
%     else
%         imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1);
%     end
% 
%     % Close the figure
%    % close(gcf);
% end





%% -----Reconstruct Vspdmd data using spDMD
% Define the time vector and the number of frames
% t = linspace(0, 42,42);
% numFrames = 42;
% 
% % Set up the figure for the animation
% figure;
% axis tight;
% xlabel('y'); % Adjust based on your data
% ylabel('z'); % Adjust based on your data
% 
% % Create a GIF file
% gifFilename = 'VspDMD_reconstruction.gif';
% %gifFilename = 'HumidityspDMD_reconstruction.gif';
% 
% for k = 1:numFrames
%     % Extract the k-th frame from Vspdmd and compute the magnitude
%     frameData = abs(reshape(Vspdmd(:, k), [40, 97])'); % Adjust dimensions if needed
%       %frameData = abs(reshape(Vspdmd(:, k)', [97, 40])); % Adjust dimensions if neede
% 
%     % Plot the current frame
%     imagesc(data.y, data.z, frameData); % Adjust axis if needed
%     colormap(brighten(redblueTecplot(21), -0.55));
%     colorbar;
%     axis xy;
%     xlim([0 2e4]);
%     ylim([0 2e4]);
%      xlabel("y");
%     ylabel("z"); 
%     title(['Time Step: ', num2str(k)]);
% 
%     % Capture the current frame
%     frame = getframe(gcf);
%     im = frame2im(frame);
%     [imind, cm] = rgb2ind(im, 256);
% 
%     % Write to GIF file
%     if k == 1
%         imwrite(imind, cm, gifFilename, 'gif', 'LoopCount', inf, 'DelayTime', 1);
%     else
%         imwrite(imind, cm, gifFilename, 'gif', 'WriteMode', 'append', 'DelayTime', 1);
%     end
% end


%% Specific the loation i-the element related to 
% x_k,i =sum_i^\infty*(\phi_j,i)*\labda_j*b_j

% Define the row and column of the point you want to highlight
% highlight_row = 31;
% highlight_col = 40;
% 
% % Calculate the linear index corresponding to this location
% highlight_index = (highlight_col - 1) * 40 + highlight_row;
% 
% for i = 1:2:rr
%     figure; % Create a new figure for each mode
%     mode_i = abs(DMDModes_xsp(:, i)); 
%     mode_i = mode_i(1:40*97)' % Extract the first 40*97 elements
% 
%     % Plot the spatial structure
%    % imagesc(data.y, data.z, reshape(mode_i, [40, 97])'); % Use transpose
%    % for correct orientation % vorticity
%     imagesc(data.y, data.z, reshape(mode_i, [97, 40])'); % humidity
%     colorbar; % Add colorbar for reference
%     colormap(brighten(redblueTecplot(21),-0.55));
%     % axis xy; % Ensure the correct orientation
%    % xlim([0 2e4]);
%    % ylim([0 2e4]);
%     xlabel("y");
%     ylabel("z"); 
%     title(['Spatial Structure of Mode ', num2str(i)]);
% 
%     % Highlight the specific point
%     hold on;
%     % Calculate the actual coordinates for y and z
%     highlight_y = data.y(highlight_row);
%     highlight_z = data.z(highlight_col);
%     plot(highlight_y, highlight_z, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'yellow'); % Use a yellow circle
%     hold off;
% end


%%-------------Subplots specific loaction

% % Define modes to plot (example: 1, 3, 5, 7, 9, 11)
% modesToPlot = 1:2:rr; % Adjust this based on your actual mode range
% numModes = length(modesToPlot);
% 
% % Number of rows and columns in the subplot grid
% numRows = 2;
% numCols = 3;
% 
% % Define the row and column of the point you want to highlight
% highlight_row = 31;
% highlight_col = 40;
% 
% % Calculate the linear index corresponding to this location
% highlight_index = (highlight_col - 1) * 40 + highlight_row;
% 
% % Create a figure for the subplots
% figure;
% 
% for k = 1:numModes
%     i = modesToPlot(k);
% 
%     % Determine subplot position
%     subplotRow = ceil(k / numCols);
%     subplotCol = mod(k - 1, numCols) + 1;
% 
%     % Extract and reshape mode_i
%     mode_i = abs(DMDModes_xsp(:, i)); 
%     mode_i = mode_i(1:40*97); % Extract the first 40*97 elements
% 
%     % Plot the spatial structure in the subplot
%     subplot(numRows, numCols, (subplotRow - 1) * numCols + subplotCol);
%     imagesc(data.y, data.z, reshape(mode_i, [40, 97])'); % Use transpose for correct orientation
%     colorbar; % Add colorbar for reference
%     colormap(brighten(redblueTecplot(21),-0.55));
%     axis xy; % Ensure the correct orientation
%     xlim([0 2e4]);
%     ylim([0 2e4]);
%     xlabel("y");
%     ylabel("z"); 
%     title(['Spatial Mode ', num2str(i)]);
% 
%     % Highlight the specific point
%     hold on;
%     % Calculate the actual coordinates for y and z
%     highlight_y = data.y(highlight_row);
%     highlight_z = data.z(highlight_col);
%     plot(highlight_y, highlight_z, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'yellow'); % Use a yellow circle
%     hold off;
% end



%%------temporal evolution

% Step 1: Define the time vector t
% t = linspace(0, 42, size(V0, 2)); % Assuming t ranges from 0 to 42 over the number of snapshots
% 
% % Step 2: Compute the time dynamics for each time step
% time_dynamics = zeros(r, length(t));
% for iter = 1:length(t)
%     time_dynamics(:, iter) = Norm_xsp .* exp(log(DEv_xsp) * t(iter));
% end
% 
% %%% DMDModes_xsp is the DMD modes (i.e., Phi) with sparsity case
% Vspdmd = DMDModes_xsp * time_dynamics;
% 
% % % Step 3: Extract the temporal evolution of the highlighted index from the original data
%  x_k_highlight = V0(highlight_index, :); % Extracting the i-th element across all snapshots
% % 
% % % Step 4: Plot the temporal evolution of the original data at the highlighted index
%  figure;
% plot(t, real(x_k_highlight), 'k', 'LineWidth', 2); % Plot the real part of the original data
% hold on; % Retain the plot for adding temporal dynamics
% 
% % Step 5: Plot the corresponding temporal dynamics from the DMD modes
% for i = 1:2:rr
%     plot(t, real(time_dynamics(i, :))); % Plot the real part of temporal dynamics for each mode
% end
% 
% % Step 6: Customize the plot
% title('Temporal Evolution at Highlighted Location and Temporal Dynamics of Modes');
% xlabel('Time');
% ylabel('Amplitude');
% legend_entries = [{'Original Data'}, arrayfun(@(x) ['Temporal Mode ', num2str(x)], 1:2:rr, 'UniformOutput', false)];
% legend(legend_entries, 'Location', 'Best'); % Add a legend to distinguish between original data and modes
% hold off;



%%%%%%%%%%%%%%%  New plot Temporal -(original, reconstruct, the decompose of ith elements)--------------------
% Define modes to plot (example: 1, 3, 5, 7, 9, 11)
% modesToPlot = 1:2:rr; % Adjust this based on your actual mode range
% numModes = length(modesToPlot);
% 
% % Number of rows and columns in the subplot grid
% numRows = 2;
% numCols = 3;
% 
% % Define the row and column of the point you want to highlight
% highlight_row = 25;
% highlight_col = 13;
% 
% % Calculate the linear index corresponding to this location
% highlight_index = (highlight_col - 1) * 40 + highlight_row;
% 
% % Step 1: Define the time vector t
% t = linspace(0, 42, size(V0, 2)); % Assuming t ranges from 0 to 42 over the number of snapshots
% 
% % Step 2: Compute the time dynamics for each time step
% time_dynamics = zeros(r, length(t));
% for iter = 1:length(t)
%     time_dynamics(:, iter) = Norm_xsp .* exp(log(DEv_xsp) * t(iter));
% end
% 
% %%% DMDModes_xsp is the DMD modes (i.e., Phi) with sparsity case
% Vspdmd = DMDModes_xsp * time_dynamics;
% 
% % Step 3: Extract the temporal evolution of the highlighted index from the original data
% %x_k_highlight = V0(highlight_index, :); % Extracting the i-th element across all snapshots
% 
% % Step 4: Extract the temporal evolution of the highlighted index from the DMD reconstruction
% Vspdmd_highlight = Vspdmd(highlight_index, :); % Extracting the i-th element across all snapshots in DMD reconstruction
% 
% % Step 5: Plot the temporal evolution of the original data at the highlighted index
% % figure;
%  plot(t,real(x_k_highlight), 'k', 'LineWidth', 2); % Plot the real part of the original data
% hold on; % Retain the plot for adding DMD reconstructed dynamics
% 
% % Step 6: Plot the reconstructed temporal evolution from DMD at the highlighted index
% plot(t, real(Vspdmd_highlight), 'r--', 'LineWidth', 2); % Plot the real part of the DMD reconstructed data
% 
% % Step 7: Plot the contributions of selected DMD modes to the temporal evolution at the highlighted index
% for i = 1:2:rr
%     % Compute the contribution of mode i to the highlight_index
%     mode_contribution = DMDModes_xsp(highlight_index, i) * time_dynamics(i, :);
% 
%     % Plot the contribution from the selected mode
%     plot(t, real(mode_contribution), '--', 'LineWidth', 1.5); % Dashed line for mode contribution
% end
% 
% % Step 8: Customize the plot
% title('Temporal Evolution at Highlighted Location: Original, DMD Reconstruction, and Mode Contributions');
% xlabel('Time');
% ylabel('Amplitude');
% legend_entries = [{'Original Data', 'DMD Reconstruction'}, arrayfun(@(x) ['Mode ', num2str(x)], 1:2:rr, 'UniformOutput', false)];
% legend(legend_entries, 'Location', 'Best'); % Add a legend to distinguish between original data, reconstruction, and modes
% hold off;
% 
% 
% 
% %===============New Temporal the i-th componet x_k,i is composed of 6
% % step 4 is cocnsist of the step 6
% modesToPlot = 1:2:rr; % Adjust this based on your actual mode range
% numModes = length(modesToPlot);
% 
% % Number of rows and columns in the subplot grid
% numRows = 2;
% numCols = 3;
% 
% % Define the row and column of the point you want to highlight
% highlight_row = 25;
% highlight_col = 13;
% 
% % Calculate the linear index corresponding to this location
% highlight_index = (highlight_col - 1) * 40 + highlight_row;
% 
% % Step 1: Define the time vector t
% t = linspace(0, 42, size(V0, 2)); % Assuming t ranges from 0 to 42 over the number of snapshots
% 
% % Step 2: Compute the time dynamics for each time step
% time_dynamics = zeros(r, length(t));
% for iter = 1:length(t)
%     time_dynamics(:, iter) = Norm_xsp .* exp(log(DEv_xsp) * t(iter));
% end
% 
% %%% DMDModes_xsp is the DMD modes (i.e., Phi) with sparsity case
% Vspdmd = DMDModes_xsp * time_dynamics;
% 
% % Step 4: Extract the temporal evolution of the highlighted index from the DMD reconstruction
% Vspdmd_highlight = Vspdmd(highlight_index, :); % Extracting the i-th element across all snapshots in DMD reconstruction
% 
% % Step 5: Plot the reconstructed temporal evolution from DMD at the highlighted index
% figure;
% plot(t, real(Vspdmd_highlight), 'r--', 'LineWidth', 2); % Plot the real part of the DMD reconstructed data
% hold on; % Retain the plot for adding DMD mode contributions
% 
% % Step 6: Plot the contributions of selected DMD modes to the temporal evolution at the highlighted index
% for i = 1:2:rr
%     % Compute the contribution of mode i to the highlight_index
%     mode_contribution = DMDModes_xsp(highlight_index, i) * time_dynamics(i, :);
% 
%     % Plot the contribution from the selected mode
%     plot(t, real(mode_contribution), '--', 'LineWidth', 1.5); % Dashed line for mode contribution
% end
% 
% % Step 7: Customize the plot
% title('Temporal Evolution at Highlighted Location: DMD Reconstruction and Mode Contributions');
% xlabel('Time');
% ylabel('Amplitude');
% legend_entries = arrayfun(@(x) ['Mode ', num2str(x)], 1:2:rr, 'UniformOutput', false);
% legend(['DMD Reconstruction', legend_entries], 'Location', 'Best'); % Add a legend to distinguish between reconstruction and modes
% hold off;
