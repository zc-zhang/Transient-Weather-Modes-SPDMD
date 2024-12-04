% Example code to compute and visualize vorticity flow from frameData
% vortex dynamics of wind field using tmpV and tmpW data
% vorflow.m
%===================
% Example code to compute and visualize vorticity flow from frameData
% vortex dynamics of wind field using tmpV and tmpW data
% vorflow.m
%===================
fig = figure;
Vorframes(121) = struct('cdata', [], 'colormap', []);

figfile = 'vorticity_flow_241202.gif';

  Vfull12=[];

for i = 1:121
    scale = 10e1; % Adjust the length of arrows for velocity components
    tmpV = squeeze(data.V(:,:,:,i)) * scale; % velocity comp in y directions  (x vertical)
    tmpW = squeeze(data.W(:,:,:,i)) * scale; % velocity comp in z directions ( y horizational)
    tmpPRES=squeeze(data.PRES(:,:,:,i)); % pressure
    tmpT=squeeze(data.T(:,:,:,i)); % temperarure

    % Ensure data.y and data.z are vectors with the correct sizes
    y = data.y(:); % 40x1 vector    % y=data.y
    z = data.z(:); % 97x1 vector    % z= data.z

    % Initialize gradients
    dW_dy = zeros(size(tmpW));
    dW_dz = zeros(size(tmpW));
    dV_dy = zeros(size(tmpV));
    dV_dz = zeros(size(tmpV));

    % Calculate dy and dz (assuming uniform spacing)
    dy = mean(diff(y));
    dz = mean(diff(z));

    % Calculate gradients using finite differences
    for j = 2:length(y)-1
        for k = 2:length(z)-1
            dW_dy(j,k) = (tmpW(j+1,k) - tmpW(j-1,k)) / (2*dy);
            dW_dz(j,k) = (tmpW(j,k+1) - tmpW(j,k-1)) / (2*dz);
            dV_dy(j,k) = (tmpV(j+1,k) - tmpV(j-1,k)) / (2*dy);
            dV_dz(j,k) = (tmpV(j,k+1) - tmpV(j,k-1)) / (2*dz);
        end
    end

    % Handle boundaries with forward/backward differences
    for j = 1
        for k = 1:length(z)
            dW_dy(j,k) = (tmpW(j+1,k) - tmpW(j,k)) / dy;
            dV_dy(j,k) = (tmpV(j+1,k) - tmpV(j,k)) / dy;
        end
    end
    for j = length(y)
        for k = 1:length(z)
            dW_dy(j,k) = (tmpW(j,k) - tmpW(j-1,k)) / dy;
            dV_dy(j,k) = (tmpV(j,k) - tmpV(j-1,k)) / dy;
        end
    end
    for j = 1:length(y)
        for k = 1
            dW_dz(j,k) = (tmpW(j,k+1) - tmpW(j,k)) / dz;
            dV_dz(j,k) = (tmpV(j,k+1) - tmpV(j,k)) / dz;
        end
    end
    for j = 1:length(y)
        for k = length(z)
            dW_dz(j,k) = (tmpW(j,k) - tmpW(j,k-1)) / dz;
            dV_dz(j,k) = (tmpV(j,k) - tmpV(j,k-1)) / dz;
        end
    end

    % Calculate vorticity components  
    % vorticity in the y-z plane (orthogonal to) x-direction
    vorticity_x = dW_dy - dV_dz; % old: % vorticity_y = dW_dy - dV_dz; 
   
   % vorticity_z = dV_dy - dW_dz; % no need to compute this term


    % % Options-1
% Calculate the stream function
    %psi = cumsum(vorticity_y, 1) * dz - cumsum(vorticity_z, 2) * dy;

     % magnitude_vorticity = abs(vorticity_x);
       magnitude_vorticity = sqrt(vorticity_x.^2);
    % Plot vorticity using contour plot
     contourf(y, z, magnitude_vorticity', 50, 'LineColor', 'none');
     % magnitude_vorticity = sqrt(vorticity_x.^2 + vorticity_z.^2);
     colorbar;
     %colormap('jet');    % intensity (or magnitude) of vorticity, with colors indicating how strong the rotational effects are.
      colormap(brighten(redblueTecplot(21),-0.55));
     hold on;

   %  Plot vorticity using quiver plot  % direction and strength of the rotational vectors,
   %  scale_vorticity = 5 * 10e1; % Adjust as necessary
     %  quiver(y, z, vorticity_x' * scale_vorticity, zeros(size(vorticity_x))', 'k');
     % %  %  quiver(y, z, vorticity_y' * scale_vorticity, vorticity_z' * scale_vorticity, 'k'); 
    % Plot streamlines to visualize rotation
      %[Y, Z] = meshgrid(y, z);
     %  startx = linspace(min(y), max(y), 20);
     %  starty = linspace(min(z), max(z), 20);
     %  streamline(Y, Z,vorticity_y', vorticity_z', startx, starty);

    hold off;


    % Customize plot 
    xlim([0 2e4]);
    ylim([0 2e4]);
    xlabel("y");
    ylabel("z");
    title(['Magnitude of Vortex at time ' num2str(i)]);
    drawnow;
    Vorframes(i) = getframe(fig);

    % Save vorticity frames for animation
    [Avor, mapvor] = rgb2ind(frame2im(Vorframes(i)), 256);
    if i == 1
        imwrite(Avor, mapvor, figfile, 'gif', 'DelayTime', 1);
    else
        imwrite(Avor, mapvor, figfile, 'gif', 'DelayTime', 1, 'WriteMode', 'append');
    end

    %    % Vectorize vorticity componentsVorframes
    % vort_y_vec = vorticity_y(:);
    % vort_z_vec = vorticity_z(:);
    % 
    % % Concatenate into a single vector for each snapshot
    % snapshot_data = [vort_y_vec; vort_z_vec];

        % Vectorize vorticity magnitude
    magnitude_vorticity_vec = magnitude_vorticity(:);

    % Append to Vfull (only the saclar field of vortex)

 Vfull12 = [Vfull12, magnitude_vorticity_vec];

  % VFull(:, i) = magnitude_vorticity(:);
end

% Save Voriticity flow data for later use in DMD
save('Vfull12.mat', 'Vfull12');


% Extract V0 and V1 from Xfull
% V0 = Vfull(:, 1:end-1);
% V1 = Vfull(:, 2:end);

% perform DMD 
% use spdmd=vorscale.m