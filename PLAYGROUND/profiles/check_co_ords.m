clear variables
close all
CameraNo = 1;
base = { ...
    'D:\Full\Processed_PIV_validation\90degree_400light_100hz_3000dt', ...
    'D:\Full\Processed_PIV_validation\90degree_250light_250hz_1000dt', ...
    };

imageCount = 18000;
caseImages = 3600;
endpoint = '';
i = 5;
plot_type = 'v';

figure;
for run = 1:length(base)
    % Load statistics
    statsloc = fullfile(base{run}, 'Statistics', num2str(imageCount), ...
                ['Cam', num2str(CameraNo)], 'Instantaneous', 'Calibrated');
    StatsData = load(fullfile(statsloc, 'MeanStats16x16.mat'));
    mean_U = StatsData.mean_U;   % Expected to be vectorized for the region
    
    % Load coordinates
    directory_path = fullfile(base{run}, 'CalibratedPIV', num2str(imageCount), ['Cam', num2str(CameraNo)], 'Instantaneous');
    Co_ords = load(fullfile(directory_path, "Co_ords"));
    
    % Extract relevant coordinates
    ycorners = [Co_ords.Co_ords(i).y(1,1), Co_ords.Co_ords(i).y(end, end)]; %, xcoords(m, 1), xcoords(m, n)];
    xcorners = [Co_ords.Co_ords(i).x(1,1), Co_ords.Co_ords(i).x(end, end)];
    
    % Create subplot
    subplot(1, length(base), run);
    imagesc(xcorners, ycorners, mean_U);  % Plot the mean velocity data
    hold on;
    
    hold off;
    axis equal;
    colorbar;
    title([run, num2str(run)]);
    xlabel('X');
    ylabel('Y');
    daspect([1 1 1])
    set(gca, 'YDir', 'normal', 'FontSize', 16, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');
end
