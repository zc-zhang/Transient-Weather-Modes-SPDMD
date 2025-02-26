% Load data
% ncfile = "your_data.nc"; % Replace with your actual file
% data.PREC = ncread(ncfile, "PREC");  % Size: 1x40x121
% data.y = ncread(ncfile, "y");        % Size: 40x1
% data.time = ncread(ncfile, "time");  % Size: 121x1

% Set up figure
figure;
%colormap(jet);
 colorbar; % Add colorbar for reference
 colormap(brighten(redblueTecplot(21),-0.05));
caxis([min(data.PREC(:)), max(data.PREC(:))]); % Adjust color limits

% Create animation
% for t = 1:121
%     imagesc(data.y, [], squeeze(data.PREC(1,:,t))); % Plot precipitation as a heatmap
%     colorbar;
%     xlabel('y (horizontal)');
%     ylabel('Precipitation Intensity');
%     title(sprintf('Precipitation Evolution at t = %.1f min', data.time(t)/60));
%     pause(0.1); % Adjust for animation speed
% end

%% Line Plot Animation (plot)
% If precipitation is 1D (along y), a simple line plot over time might be clearer:
figure;
for t = 1:121
    plot(data.y, squeeze(data.PREC(1,:,t)), 'LineWidth', 2);
    ylim([min(data.PREC(:)), max(data.PREC(:))]); % Keep consistent scale
    xlabel('y (horizontal)');
    ylabel('Precipitation Intensity');
    title(sprintf('Precipitation Evolution at t = %.1f min', data.time(t)/60));
    grid on;
    pause(0.1);
end



%  Bar Chart (bar)
% if you want a discrete representation of precipitation changes:
figure;
for t = 1:121
    bar(data.y, squeeze(data.PREC(1,:,t)));
    ylim([min(data.PREC(:)), max(data.PREC(:))]); % Keep consistent scale
    xlabel('y (horizontal)');
    ylabel('Precipitation Intensity');
    title(sprintf('Precipitation at t = %.1f min', data.time(t)/60));
    pause(0.1);
end



%% Waterfall Plot (waterfall) 
% This shows precipitation evolution over time in 3D:
figure;
waterfall(data.y, data.time/60, squeeze(data.PREC(1,:,:))');
xlabel('y (horizontal)');
ylabel('Time (minutes)');
zlabel('Precipitation Intensity');
title('Precipitation Evolution');
%colormap(jet);
colormap(brighten(redblueTecplot(21),-0.45));
colorbar;

%% Surface Plot (surf)
% This can give a smoother look over time:

figure;
surf(data.y, data.time/60, squeeze(data.PREC(1,:,:))');
shading interp;
xlabel('y (horizontal)');
ylabel('Time (minutes)');
zlabel('Precipitation Intensity');
title('Precipitation Evolution');
%colormap(jet);
colormap(brighten(redblueTecplot(21),-0.45));
colorbar;
view(2); % Top-down for better readability


