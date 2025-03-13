% Initialize figure
fig = figure;
% Define GIF filename for humidity ratio
humidity_giffile = 'humidity_evolution_250313.gif';

% Loop through time steps
for i = 1:length(data.time)
    % Extract humidity ratio data for current time step
    tmpQV = squeeze(data.QV(:,:,:,i)); % Assuming data.QV is your humidity data
    
    % Plot humidity ratio
    imagesc(flip(tmpQV'));
    colorbar;
    colormap(brighten(redblueTecplot(21), -0.05));
    xlabel('y');
    ylabel('z');
   % title(['Humidity Ratio (g/kg) at t = ' num2str(data.time(i)/60) ' min']);
   title(['Humidity Ratio (g/kg) t= ' num2str(i) ]);
    drawnow;
    
    % Capture the frame
    frame = getframe(fig);
    [A, map] = rgb2ind(frame2im(frame), 256);
    
    % Write to GIF file
    if i == 1
        imwrite(A, map, humidity_giffile, 'gif', 'DelayTime', 0.1, 'Loopcount', 1);
    else
        imwrite(A, map, humidity_giffile, 'gif', 'DelayTime', 0.1, 'WriteMode', 'append');
    end
end

% Close the figure (optional)
close(fig);