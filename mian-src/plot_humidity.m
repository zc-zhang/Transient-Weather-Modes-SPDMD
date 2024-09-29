% plot for humidity data

% sort amplitudes for xpol, xsp, xdmd

[Norm_xdmd,Index_xdmd] = sort(abs(xdmd),'descend');
DEv_xdmd = Edmd(Index_xdmd);   %discrete-eigenvalues 
DMDModes_xdmd = Phi(:,Index_xdmd);

kk=110;
rr=answer.Nz(kk);

%sort_xsp = answer.xsp(:,kk)  
% sort of the large amplitudes rather than rank of Eigvals
[Norm_xsp,Index_xsp] = sort(abs(answer.xsp(:,kk)),'descend');  %return the value of order/peak ofamplitudes
DEv_xsp = Edmd(Index_xsp);   %discrete-eigenvalues 
DMDModes_xsp = Phi(:,Index_xsp);


[Norm_xpol,Index_xpol] = sort(abs(answer.xpol(:,kk)),'descend');
DEv_xpol = Edmd(Index_xpol);   %discrete-eigenvalues 
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
     colormap(brighten(redblueTecplot(21),-0.55));
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

% Define the row and column of the point you want to highlight, in[40,97]
 % highlight_row = 30;  %  <40
 % highlight_col = 22;   % <97

  highlight_col = 30;  %  <40
  highlight_row = 25;   % <97

% Calculate the linear index corresponding to this location
%highlight_index = (highlight_col - 1) * 40 + highlight_row;
highlight_index = (highlight_col - 1) * 97 + highlight_row;
%highlight_index = (highlight_row - 1) * 97 + highlight_col;
% Create a figure for the subplots
figure;

for k = 1:numModes
    i = modesToPlot(k);
    
    % Determine subplot position
    subplotRow = ceil(k / numCols);
    subplotCol = mod(k - 1, numCols) + 1;
    
    % Extract and reshape mode_i
    %mode_i = abs(DMDModes_xdmd(:, i)); % DMD modes
    %mode_i = abs(DMDModes_xsp(:, i)); spDMDModes
    %mode_i = abs(CKModes(:, i));
    mode_i = mode_i(1:40*97); % Extract the first 40*97 elements
    
    % Plot the spatial structure in the subplot
    subplot(numRows, numCols, (subplotRow - 1) * numCols + subplotCol);
    imagesc(data.y, data.z, reshape(mode_i, [97, 40])); % Use transpose for correct orientatio  
    colorbar; % Add colorbar for reference
    colormap(brighten(redblueTecplot(21),-0.55));
    xlim([0 2e4]);
    ylim([0 2e4]);
    xlabel("y");
    ylabel("z"); 
    title(['Spatial Mode ', num2str(i)]);
    
    % Highlight the specific point
    hold on;
    % Calculate the actual coordinates for y and z
    % highlight_y = data.y(highlight_row);
    % highlight_z = data.z(highlight_col);
    highlight_y = data.y(highlight_col);
    highlight_z = data.z(highlight_row);
    plot(highlight_y, highlight_z, 'ko', 'MarkerSize', 9, 'MarkerFaceColor', 'yellow'); % Use a yellow circle
    hold off;
end


% sort or account the number of zero amplitude related eigenvals 
 nonzero_indices = find(answer.xsp(:,kk)); % Step 1: Find nonzero amplitudes
  [~, sort_order] = sort(nonzero_indices, 'descend'); 
  ival = nonzero_indices(sort_order);    % also ralted to `DEv_xsp' (the order of eigvals)

 
% plot eigenvalues vs spdmd
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

% t = linspace(0, 39,39);
% numFrames = 40;
% 
% % Set up the figure for the animation
% figure;
% axis tight;
% xlabel('y'); % Adjust based on your data
% ylabel('z'); % Adjust based on your data
% 
% % Create a GIF file
% %gifFilename = 'VspDMD_reconstruction.gif';
% gifFilename = 'HumidityspDMD_reconstruction.gif';
% 
% for k = 1:numFrames
%     % Extract the k-th frame from Vspdmd and compute the magnitude
%       frameData = abs(reshape(Vspdmd(:, k)', [97, 40])); % humidity
% 
%     % Plot the current frame
%     imagesc(data.y, data.z, frameData); % Adjust axis if needed
%     colormap(brighten(redblueTecplot(21), -0.55));
%     colorbar;
%     %axis xy;
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