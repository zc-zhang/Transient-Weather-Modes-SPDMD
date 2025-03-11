%% comp_scalar_flow.m

% generate the magnitude of scalar filed of data.z and data.y

%%%%%%%%%%%% generate the Snapshot data matrices %==================



fig = figure;
frames(121) = struct('cdata', [], 'colormap', []);

figfile = 'scalar_counter_imag_250312.gif';

% Define spatial dimensions
ny = length(data.y); % Number of y-points
nx = length(data.z); % Number of x-points
M2 = 121;  % Number of snapshots (adjust as needed)

% Initialize the snapshot matrix
%XMaFull = [];


for i = 1:M2
    scale = 1e0; % Adjust the length of the arrows
    tmpV = squeeze(data.V(:,:,:,i)) * scale;
    tmpW = squeeze(data.W(:,:,:,i)) * scale;

    % Calculate the magnitude of the vector field
    magnitude = sqrt(tmpV.^2 + tmpW.^2);
    
    % Plot magnitude using imagesc
  % imagesc(data.y, data.z, magnitude')

   contourf(data.y, data.z, sqrt(magnitude)', 50, 'LineColor', 'none');
   
     colormap(brighten(redblueTecplot(21),-0.45))
     %colormap(jet); % Choose a colormap
    % colorbar; % Add colorbar for reference
    axis xy; % Ensure the correct orientation
    xlim([0 2e4]);
    ylim([0 2e4]);
    xlabel("y");
    ylabel("z"); 
    title(['Magnitude of scalar field at time ' num2str(i)]);
    drawnow;
    
    % Capture the frame for GIF
    frames(i) = getframe(fig);
    [A, map] = rgb2ind(frame2im(frames(i)), 256);
    if i == 1
        imwrite(A, map, figfile, 'gif', 'DelayTime', 1);
    else
        imwrite(A, map, figfile, 'gif', 'DelayTime', 1, 'WriteMode', 'append');
    end

    % Calculate the magnitude of scalar 
   % magnitude = sqrt(tmpV.^2 + tmpW.^2);

    % Flatten the 2D magnitude matrix into a 1D vector
   % XMaFull(:, i) = magnitude(:);

end

% Save Scalar flow data for later use in DMD
%save('XMaFull.mat', 'XMaFull');
