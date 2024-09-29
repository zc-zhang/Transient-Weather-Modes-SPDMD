%
% datahandle0.m -  by YS 20240501
%
%ncfile = "./ideal/merged_history.pe000000.nc";

ncfile="D:\Susuki Lab\Testing_Code\data-weather\Ensemble SCALE data Test\data20240805b\ideal/merged_history.pe000000.nc";

data=struct('time','x','y','z','V','W','PREC',[]); %'QV','PRES','T',[]);

data.time=ncread(ncfile,"time");
data.x=ncread(ncfile,"x");
data.y=ncread(ncfile,"y");
data.z=ncread(ncfile,"z");
data.V=ncread(ncfile,"V"); % Y方向水平速度成分(3次元) 
data.W=ncread(ncfile,"W"); % 鉛直方向速度成分(3次元)
data.PREC=ncread(ncfile,"PREC"); % 降水強度(2次元)
data.QV=ncread(ncfile,"QV"); % 水蒸気比湿  specific humidity or humidity ratio
data.PRES=ncread(ncfile,"PRES"); % 圧力 % pressure
data.T=ncread(ncfile,"T"); % 温度 % temperature

save("scaledata240805b.mat","data")

fig=figure;
frames(121) = struct('cdata', [], 'colormap', [])

%figfile = 'windfield240805b.gif';
figfile = 'windfield240928c.gif';

Pdata =[]; 
P_rms =[]; % pressure root mean square (RMS) error over time
Tdata=[]; 
Hdata=[];


for i=1:1:length(data.time);
    scale=10e1; % 矢印の長さの調整
    tmpV=squeeze(data.V(:,:,:,i))*scale;
    tmpW=squeeze(data.W(:,:,:,i))*scale;
    tmpPRES=squeeze(data.PRES(:,:,:,i));
    tmpT=squeeze(data.T(:,:,:,i));
    % optional:
    tmpQV=squeeze(data.QV(:,:,:,i));
    %tmpPREC=squeeze(data.PREC(:,:,:,i));
   

    subplot(2,2,1);
    h=quiver(data.y,data.z,tmpV',tmpW','off');
    xlim([0 2e4]);
    ylim([0 2e4]);
    xlabel("y");
    ylabel("z");
    title(['time(in sec) ' num2str(data.time(i))]);
    drawnow;
    frames(i)=getframe(fig);
    
    [A,map]=rgb2ind(frame2im(frames(i)),256);
    if i==1
        imwrite(A,map,figfile,'gif','DelayTime',1);
    else
        imwrite(A,map,figfile,'gif','DelayTime',1,'WriteMode','append');
    end

    subplot(2,2,2);
    imagesc(flip(tmpPRES'));
    colorbar;
    title(['pressure (in Pa)']);
    drawnow;
   

    subplot(2,2,3);
    imagesc(flip(tmpT'));
    colorbar;
    title(['temperature (in K)']);
    drawnow;

    subplot(2,2,4);
    imagesc(flip(tmpQV'));
    colorbar;
    title(['humidity ratio (in g/kg)']); % or specific humidity
    drawnow;
 
    % vectorize the data (stcak data matricies as vector)
    flip_Pdata_i=flip(tmpPRES');
    flip_Tdata_i=flip(tmpPRES');
    flip_Hdata_i=flip(tmpQV');
    
    Pdata = [Pdata,flip_Pdata_i(:)];
    Tdata = [Tdata, flip_Tdata_i(:)];
    Hdata = [Hdata, flip_Hdata_i(:)];
end

save('Pdata.mat','Pdata');   % pressure data
 save('Tdata.mat','Tdata');  % temperature data
  save('Hdata.mat','Hdata');  % temperature data

  % plot the P_RMS root mean squre 
% File name for the GIF
% gif_filename = 'RMS_dynamic_pressure.gif';
% 
% % Initialize variables
% nRows = 40;  % Number of rows in grid
% nCols = 97;  % Number of columns in grid
% num_time_steps = 121;  % Number of time steps
% 
% % Calculate the mean pressure for each grid point across all time steps
% P_mean = mean(Pdata, 2);  % Mean across time (columns), resulting in a 3880x1 vector
% 
% % Initialize matrix for storing RMS error over all time steps
% RMS_error = zeros(3880, num_time_steps);
% 
% % Calculate RMS error at each grid point for each time step
% for t = 1:num_time_steps
%     Pdata_t = Pdata(:, t);  % Pressure data at time step t (3880x1 vector)
% 
%     % Calculate RMS error for each grid point at this time step
%     RMS_error(:, t) = sqrt((Pdata_t - P_mean).^2);  
% end
% 
% % Loop over time steps to generate frames for the GIF
% for t = 1:num_time_steps
%     % Reshape the RMS error data for the chosen time step into a 40x97 grid
%     RMS_error_grid = reshape(RMS_error(:, t), nRows, nCols);
% 
%     % Plot the spatial distribution of RMS error
%     imagesc(flip(RMS_error_grid'));  % Flip and plot the grid
%     colorbar;
%     title(['RMS Error of Dynamic Pressure at Time Step ', num2str(t)]);
%     xlabel('X');
%     ylabel('Y');
% 
%     % Capture the plot as a frame
%     frame = getframe(gcf);
%     im = frame2im(frame);
%     [imind, cm] = rgb2ind(im, 256);  % Convert to indexed image
% 
%     % Write to the GIF file
%     if t == 1
%         imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
%     else
%         imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
%     end
% 
%     pause(0.1);  % Optional: Pause to control the animation speed in real time
% end
