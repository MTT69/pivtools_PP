close all;
clear variables;

%% Parameters and Base Directories
N = 18000;             % Total number of data points per directory
idx_start = 1;
idx_end = N;           % Adjust these indices if you want a different segment

% List of directories to process
base_dirs = { ...
    'E:\Full\Processed_PIV_validation\90degree_400light_100hz_3000dt', ...
    'E:\Full\Processed_PIV_validation\90degree_250light_250hz_1000dt', ...
    'E:\Full\Processed_PIV_validation\60degree_400light_100hz_3000dt_reverse', ...
    'E:\Full\Processed_PIV_validation\60degree_250light_250hz_1000dt_reverse', ...
    'E:\Full\Processed_PIV_validation\60degree_400light_100hz_3000dt', ...
    'E:\Full\Processed_PIV_validation\60degree_250light_250hz_1000dt', ...
    'E:\Full\Processed_PIV_validation\30degree_400light_100hz_3000dt_reverse', ...
    'E:\Full\Processed_PIV_validation\30degree_250light_250hz_1000dt_reverse', ...
    'E:\Full\Processed_PIV_validation\30degree_400light_100hz_3000dt', ...
    'E:\Full\Processed_PIV_validation\30degree_250light_250hz_1000dt' ...
    };

numBaseDirs = length(base_dirs);

% Define lag values (in samples): five negative steps, zero, and five positive steps.
% For a step size of 20, the lags will be: -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100.
lags = -100:20:100;
numLags = length(lags);

%% Process Each Directory
for d = 1:numBaseDirs
    base_dir = base_dirs{d};
    fprintf('Processing directory %d of %d:\n   %s\n', d, numBaseDirs, base_dir);
    
    % Construct the statistics folder and load the data
    statistics_dir = fullfile(base_dir, 'Statistics', num2str(N), 'Cam1', 'Instantaneous', 'Calibrated');
    data = load(fullfile(statistics_dir, 'cavity_stats.mat'));
    
    % Extract the flux (column 6) and delta star (column 2) signals and fill missing values
    flux6      = fillmissing(data.v_flux_all(idx_start:idx_end, 6), 'linear');
    delta_star = fillmissing(data.deltaStar_all(idx_start:idx_end, 2), 'linear');
    
    % Create one figure for the current directory
    figure;
    sgtitle(sprintf('Scatter Plots for %s', base_dir)); % Overall title for the figure
    
    % Determine subplot grid size: 4 columns per row
    ncols = 4;
    nrows = ceil(numLags/ncols);
    
    %% Loop Over Lag Values and Create Subplots
    for iLag = 1:numLags
        lag = lags(iLag);
        
        % Define overlapping indices based on lag:
        % Positive lag: flux is shifted forward relative to delta_star.
        if lag > 0
            ds = delta_star(1:end-lag);
            fl = flux6(1+lag:end);
        elseif lag < 0  % Negative lag: flux is shifted backward relative to delta_star.
            ds = delta_star(1 - lag:end);
            fl = flux6(1:end+lag);
        else  % Zero lag: use the signals as they are.
            ds = delta_star;
            fl = flux6;
        end
        
        % Create the subplot
        subplot(nrows, ncols, iLag);
        scatter(ds, fl, '.');
        title(sprintf('Lag = %d samples', lag));
        xlabel('\delta^*');
        ylabel('Flux');
    end
end
