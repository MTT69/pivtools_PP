% fft_A_P_welch needs to be run first

i = 5;
source = 'E:\Processed_PIV\90degree_400light_100hz_3000dt\CalibratedPIV\18000\Cam1\Instantaneous\RPCA';
filename = sprintf('%s\\%05d.mat', source, 1);
load(filename);
ux = piv_result(i).ux; % X-direction velocity
uy = piv_result(i).uy; % Y-direction velocity
close all
close(v);
load(fullfile(source,'Co_ords.mat'))
% Sampling rate in Hz
sampling_rate = 100;  % Data acquired at 100 Hz

t = 18000;
lower_top = -0.02; upper_top = 0.02;
lower_bottom = -0.04; upper_bottom = 0.04;


type = 'v'; % u' v'
p_t = 'peak';  % peak or trough
mp4_filename = fullfile(pwd, ['FlowDevelopment-', p_t, '-', type, '.mp4']);
try
    A;
catch
    load('E:\Processed_PIV\90degree_400light_100hz_3000dt\Statistics\18000\Cam1\Instantaneous\Calibrated\RPCA\Below\POD_stats_319x149.mat');
end
clear A_s LAM_s PHI PHI_s lambda_s ilam_s

mode = 1;
source = 'E:\Processed_PIV\90degree_400light_100hz_3000dt\CalibratedPIV\18000\Cam1\Instantaneous\RPCA';
skip = 5;  % Skip every nth row vector in making the quiver plot
step = 2; % skip every second nth column in making the quiver plot
dividing_row = 75; % vector location to split quiver plot
xlim_l = -80;
xlim_u = 80;
Frames = 41; % number of frames per peak of trough
ylim_u = 20;
ylim_l = -45;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(p_t, 'peak')
    load("peak_indices.mat")
    combined_indices = extended_peak_locations;
elseif strcmp(p_t, 'trough')
    load("trough_indices.mat")
    combined_indices = extended_trough_locations;
end

if strcmp(type, 'u')
    load("mean_U.mat");
    load("mean_V.mat")
    mean_field = mean_ux;
elseif strcmp(type, 'v')
    load("mean_V.mat");
    load("mean_U.mat");
    mean_field = mean_uy;
end


ycorners_top = [Co_ords(i).y(1, 1), Co_ords(i).y(dividing_row, end)];
xcorners_top = [Co_ords(i).x(1, 1), Co_ords(i).x(dividing_row, end)];
ycorners_bottom = [Co_ords(i).y(dividing_row, 1), Co_ords(i).y(end, end)];
xcorners_bottom = [Co_ords(i).x(dividing_row, 1), Co_ords(i).x(end, end)];



Total_flow_top = zeros(size(ux(1:dividing_row,:)));
Total_flow_bottom = zeros(size(ux(dividing_row:end,:)));  
Total_flow_top_ux = zeros(size(ux(1:step:dividing_row, 1:step:end)));
Total_flow_top_uy = Total_flow_top_ux;
Total_flow_bottom_ux = zeros(size(ux(dividing_row+1:step:end, 1:step:end)));
Total_flow_bottom_uy = Total_flow_bottom_ux;


%%
numCells = length(combined_indices) / Frames;

% This is preallocation for the shorter peak average video
developmentVideo_top = cell(1, numCells);
developmentVideo_bottom = cell(1, numCells);
developmentVideo_top_ux = cell(1, numCells);
developmentVideo_top_uy = cell(1, numCells);
developmentVideo_bottom_ux = cell(1, numCells);
developmentVideo_bottom_uy = cell(1, numCells);
% Initialize each cell with the relevant size array

for k = 1:numCells
    developmentVideo_top{k} = zeros(size(Total_flow_top)); % Match Total_flow_top dimensions
    developmentVideo_bottom{k} = zeros(size(Total_flow_bottom)); % Match Total_flow_bottom dimensions
    developmentVideo_top_ux{k} = zeros(size(ux(1:step:dividing_row, 1:step:end)));
    developmentVideo_top_uy{k} = zeros(size(ux(1:step:dividing_row, 1:step:end)));
    developmentVideo_bottom_ux{k} = zeros(size(ux(dividing_row+1:step:end, 1:step:end)));
    developmentVideo_bottom_uy{k} = zeros(size(ux(dividing_row+1:step:end, 1:step:end)));
end

%% 
% video setup

v = VideoWriter(mp4_filename, 'MPEG-4');  % Create a VideoWriter object for MP4
v.FrameRate = 3;  % Set the frame rate (e.g., 10 frames per second)
open(v);  % Open the video writer
count = 0;
index = 0;
for imNo = combined_indices
    % Load velocity field data for the current frame
    index = index+1;
    filename = sprintf('%s\\%05d.mat', source, imNo);
    load(filename);
    ux = piv_result(i).ux;
    uy = piv_result(i).uy;

    if strcmp(type, 'u')
        fluc = ux - mean_field;
    elseif strcmp(type, 'v')
        fluc = uy - mean_field;
    end

    % total raw velocity stepped
    Total_flow_top_ux = Total_flow_top_ux + ux(1:step:dividing_row, 1:step:end);
    Total_flow_top_uy = Total_flow_top_uy + uy(1:step:dividing_row, 1:step:end);
    Total_flow_bottom_ux = Total_flow_bottom_ux + ux(dividing_row+1:step:end, 1:step:end);
    Total_flow_bottom_uy = Total_flow_bottom_uy + uy(dividing_row+1:step:end, 1:step:end);

    % total fluctating 
    fluc_top = fluc(1:dividing_row,:);
    fluc_bottom = fluc(dividing_row:end,:);  
    Total_flow_top = Total_flow_top + fluc_top;
    Total_flow_bottom = Total_flow_bottom + fluc_bottom;

    % Average accross peaks for secondary video

    if mod(index - (Frames+1) ,Frames) == 0 && imNo >= Frames
        count = count +1;
    end
    developmentVideo_bottom{count} = developmentVideo_bottom{count} + fluc_bottom;
    developmentVideo_top{count} = developmentVideo_top{count} + fluc_top;
   


    b_mask = piv_result(i).b_mask;
    % Create a new figure
    fig = figure;
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);


    mask = double(b_mask);
    mask(mask==1) =inf;
    ax3 = axes;
    j = imagesc(xcorners_bottom,ycorners_bottom, mask(dividing_row:end,:));
    set(j, 'AlphaData', isinf(mask(dividing_row:end,:)));
    set(gca, 'YDir', 'normal','FontSize',18);
    colormap(ax3, [0 0 0]);
    % Align axes

    % Make sure the figure aspect ratio is consistent
    daspect([1 1 1]);   

    % Hide the second set of axes

    ax3.XTick = [];
    ax3.YTick = [];
    ax3.Visible = 'off';

    ax2 = axes; % Create second axes
    imagesc(ax2, xcorners_bottom, ycorners_bottom, fluc_bottom, [lower_bottom, upper_bottom]);
    daspect(ax2, [1 1 1]);
    set(ax2, 'YDir', 'normal', 'FontSize', 18);
    daspect([1 1 1]);
    linkaxes([ax2, ax3]);
    custommap = redbluezero(lower_bottom, upper_bottom);
    colormap(ax2, custommap);
    cbar2 = colorbar(ax2, 'eastoutside');
    cbar2.Label.String = 'Bottom Data'; % Add label to colorbar if needed


    ax3.Position = ax2.Position;  
    uistack(ax2, 'bottom'); % Push ax2 to the bottom of the stack

    % Remove extra space between plots
    ax3.XLim = ax2.XLim;
    ax3.YLim = ax2.YLim;


    ax1 = axes; % Create first axes
    imagesc(ax1, xcorners_top, ycorners_top, fluc_top, [lower_top, upper_top]);
    daspect(ax1, [1 1 1]);
    set(ax1, 'YDir', 'normal', 'FontSize', 18);
    daspect([1 1 1]);
    cbar1 = colorbar(ax1, 'westoutside');
    cbar1.Label.String = 'Top Data'; % Add label to colorbar if needed

    ax2.Position = ax1.Position;          % Match the position for overlay
    ax1.Visible = 'off';
    % Link axes for synchronized zooming and panning
    linkaxes([ax1, ax2]);
    new_position = [0.1, 0.4, 0.8, 0.55];

    % Set the position for each axes
    ax1.Position = new_position;
    ax2.Position = new_position;
    ax3.Position = new_position;

    hold on;


    % Top region: rows 1 to dividing_row
    x_top = Co_ords(i).x(1:step:dividing_row, 1:step:end);
    y_top = Co_ords(i).y(1:step:dividing_row, 1:step:end);
    ux_top = ux(1:step:dividing_row, 1:step:end);
    uy_top = uy(1:step:dividing_row, 1:step:end);
    developmentVideo_top_ux{count} = developmentVideo_top_ux{count}+ ux_top;
    developmentVideo_top_uy{count} = developmentVideo_top_uy{count}+ uy_top;


    hQuiverTop = quiver(x_top, y_top, ux_top*5, uy_top*5, "off",  'k');  % Black arrows for top region



    % Bottom region: rows dividing_row+1 to end
    x_bottom = Co_ords(i).x(dividing_row+1:step:end, 1:step:end);
    y_bottom = Co_ords(i).y(dividing_row+1:step:end, 1:step:end);
    ux_bottom = ux(dividing_row+1:step:end, 1:step:end);
    uy_bottom = uy(dividing_row+1:step:end, 1:step:end);
    developmentVideo_bottom_ux{count} = developmentVideo_bottom_ux{count}+ ux_bottom;
    developmentVideo_bottom_uy{count} = developmentVideo_bottom_uy{count}+ uy_bottom;

    hQuiverBottom = quiver(x_bottom, y_bottom, ux_bottom*40, uy_bottom*40, "off", 'r');  % Red arrows for bottom region

    xlim([xlim_l xlim_u])
    ylim([ylim_l ylim_u])


    ax4 = axes('Position', [0.1, 0.05, 0.8, 0.25]);  % Adjust position as needed

    time_coefficient = A(1:t, mode);

    % Plot the temporal coefficients
    plot((1:t) * dt, time_coefficient, 'LineWidth', 1.5);
    % Right subplot: Temporal coefficient plot with red marker

    hold on;

    % Mark the current time instance with a red dot
    current_time = imNo * dt;  % Compute the current time in seconds
    plot(current_time, time_coefficient(imNo), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
    hold off;

    % Label axes and title
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Temporal Coefficients (Current Instance Marked)');
    grid on;



    % Capture the current frame
    frame = getframe(fig);

    % Write the frame to the video
    writeVideo(v, frame);

    % Close the current figure
    close(fig);
end
% average fluctuating component
Average_flow_bottom = Total_flow_bottom/length(combined_indices); 
Average_flow_top = Total_flow_top/length(combined_indices);

% Co-ords
x_top = Co_ords(i).x(1:step:dividing_row, 1:step:end);
y_top = Co_ords(i).y(1:step:dividing_row, 1:step:end);
x_bottom = Co_ords(i).x(dividing_row+1:step:end, 1:step:end);
y_bottom = Co_ords(i).y(dividing_row+1:step:end, 1:step:end);


% average velocities
averageUx_top = Total_flow_top_ux/length(combined_indices);
averageUy_top = Total_flow_top_uy/length(combined_indices);
averageUx_bottom = Total_flow_bottom_ux/length(combined_indices);
averageUy_bottom = Total_flow_bottom_uy/length(combined_indices);


close(v);


%% create mean flow field


b_mask = piv_result(i).b_mask;
% Create a new figure
fig = figure;
set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);


mask = double(b_mask);
mask(mask==1) =inf;
ax3 = axes;
j = imagesc(xcorners_bottom,ycorners_bottom, mask(dividing_row:end,:));
set(j, 'AlphaData', isinf(mask(dividing_row:end,:)));
set(gca, 'YDir', 'normal','FontSize',18);
colormap(ax3, [0 0 0]);
% Align axes

% Make sure the figure aspect ratio is consistent
daspect([1 1 1]);   

% Hide the second set of axes

ax3.XTick = [];
ax3.YTick = [];
ax3.Visible = 'off';

ax2 = axes; % Create second axes
imagesc(ax2, xcorners_bottom, ycorners_bottom, Average_flow_bottom, [lower_bottom, upper_bottom]);
daspect(ax2, [1 1 1]);
set(ax2, 'YDir', 'normal', 'FontSize', 18);
daspect([1 1 1]);
linkaxes([ax2, ax3]);
% Set proper positioning of axes
custommap = redbluezero(lower_bottom, upper_bottom);
colormap(ax2, custommap);

ax3.Position = ax2.Position;  
uistack(ax2, 'bottom'); % Push ax2 to the bottom of the stack

% Remove extra space between plots
ax3.XLim = ax2.XLim;
ax3.YLim = ax2.YLim;


ax1 = axes; % Create first axes
imagesc(ax1, xcorners_top, ycorners_top, Average_flow_top, [lower_top, upper_top]);
daspect(ax1, [1 1 1]);
set(ax1, 'YDir', 'normal', 'FontSize', 18);
daspect([1 1 1]);

ax2.Position = ax1.Position;          % Match the position for overlay
ax1.Visible = 'off';
% Link axes for synchronized zooming and panning
linkaxes([ax1, ax2]);


hold on;


hQuiverTop = quiver(x_top, y_top, averageUx_top*5, averageUy_top*5, "off",  'k');  % Black arrows for top region


hQuiverBottom = quiver(x_bottom, y_bottom, averageUx_bottom*30, averageUy_bottom*30, "off", 'r');  % Red arrows for bottom region


daspect([1 1 1]);

xlim([xlim_l xlim_u])
ylim([ylim_l ylim_u]);


set(gca, 'YDir', 'normal', 'FontSize', 18);


% Define the file name for saving the figure
file_name = fullfile(pwd, ['AverageFlow-', p_t, '-', type]);

% Save the figure as a JPEG file
saveas(gcf, [file_name, '.jpeg']);

% Save the figure as a MATLAB figure file
savefig(gcf, [file_name, '.fig']);

%% create flow development video

mp4_filename = fullfile(pwd, ['FlowDevelopment-average', p_t, '-', type, '.mp4']);

v = VideoWriter(mp4_filename, 'MPEG-4');  % Create a VideoWriter object for MP4
v.FrameRate = 3;  % Set the frame rate (e.g., 10 frames per second)
open(v);  % Open the video writer
for k = 1:numCells
    developmentVideo_top{k} = developmentVideo_top{k} / Frames;
    developmentVideo_bottom{k} = developmentVideo_bottom{k} / Frames;
    developmentVideo_top_ux{k} = developmentVideo_top_ux{k} / Frames;
    developmentVideo_top_uy{k} = developmentVideo_top_uy{k} / Frames;
    developmentVideo_bottom_ux{k} = developmentVideo_bottom_ux{k} / Frames;
    developmentVideo_bottom_uy{k} = developmentVideo_bottom_uy{k} / Frames;
end

for L = 1:numCells

    b_mask = piv_result(i).b_mask;
    % Create a new figure
    fig = figure;
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);


    mask = double(b_mask);
    mask(mask==1) =inf;
    ax3 = axes;
    j = imagesc(xcorners_bottom,ycorners_bottom, mask(dividing_row:end,:));
    set(j, 'AlphaData', isinf(mask(dividing_row:end,:)));
    set(gca, 'YDir', 'normal','FontSize',18);
    colormap(ax3, [0 0 0]);
    % Align axes

    % Make sure the figure aspect ratio is consistent
    daspect([1 1 1]);   

    % Hide the second set of axes

    ax3.XTick = [];
    ax3.YTick = [];
    ax3.Visible = 'off';

    ax2 = axes; % Create second axes
    imagesc(ax2, xcorners_bottom, ycorners_bottom, developmentVideo_bottom{L}, [lower_bottom, upper_bottom]);
    daspect(ax2, [1 1 1]);
    set(ax2, 'YDir', 'normal', 'FontSize', 18);
    daspect([1 1 1]);
    linkaxes([ax2, ax3]);
    custommap = redbluezero(lower_bottom, upper_bottom);
    colormap(ax2, custommap);
    cbar2 = colorbar(ax2, 'eastoutside');
    cbar2.Label.String = 'Bottom Data'; % Add label to colorbar if needed


    ax3.Position = ax2.Position;  
    uistack(ax2, 'bottom'); % Push ax2 to the bottom of the stack

    % Remove extra space between plots
    ax3.XLim = ax2.XLim;
    ax3.YLim = ax2.YLim;


    ax1 = axes; % Create first axes
    imagesc(ax1, xcorners_top, ycorners_top, developmentVideo_top{L}, [lower_top, upper_top]);
    daspect(ax1, [1 1 1]);
    set(ax1, 'YDir', 'normal', 'FontSize', 18);
    daspect([1 1 1]);
    cbar1 = colorbar(ax1, 'westoutside');
    cbar1.Label.String = 'Top Data'; % Add label to colorbar if needed

    ax2.Position = ax1.Position;          % Match the position for overlay
    ax1.Visible = 'off';
    % Link axes for synchronized zooming and panning
    linkaxes([ax1, ax2]);
    new_position = [0.1, 0.4, 0.8, 0.55];

    % Set the position for each axes
    ax1.Position = new_position;
    ax2.Position = new_position;
    ax3.Position = new_position;
   
    hold on;


    % Top region: rows 1 to dividing_row
    x_top = Co_ords(i).x(1:step:dividing_row, 1:step:end);
    y_top = Co_ords(i).y(1:step:dividing_row, 1:step:end);


    hQuiverTop = quiver(x_top, y_top, developmentVideo_top_ux{L}*5, developmentVideo_top_uy{L}*5, "off",  'k');  % Black arrows for top region



    % Bottom region: rows dividing_row+1 to end
    x_bottom = Co_ords(i).x(dividing_row+1:step:end, 1:step:end);
    y_bottom = Co_ords(i).y(dividing_row+1:step:end, 1:step:end);
     

    hQuiverBottom = quiver(x_bottom, y_bottom, developmentVideo_bottom_ux{L}*40, developmentVideo_bottom_uy{L}*40, "off", 'r');  % Red arrows for bottom region

    xlim([xlim_l xlim_u])
    ylim([ylim_l ylim_u])

    % Capture the current frame
    frame = getframe(fig);

    % Write the frame to the video
    writeVideo(v, frame);

    % Close the current figure
    close(fig);
end

close(v)
