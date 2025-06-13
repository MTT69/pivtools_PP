%% averaged_crosscorr_model.m
% This script processes data from multiple base directories.
% For each directory, the 18,000 data points are split into 5 segments 
% (each of 3,600 points). For each segment the cross‐correlation is computed 
% for a set of signals and then averaged over the segments.
% Plots are produced for each base directory separately.

close all;
clear variables;

%% Parameters and Base Directories
N = 18000;             % Total number of data points per directory
segmentLength = 3600;  % Number of points per segment
nSegs = N / segmentLength;  % Number of segments per directory (should be 5)

base_dirs = { ...
%     'D:\Full\Processed_PIV_validation\90degree_400light_100hz_3000dt', ...
%     'D:\Full\Processed_PIV_validation\90degree_250light_250hz_1000dt', ...
%     'D:\Full\Processed_PIV_validation\60degree_400light_100hz_3000dt_reverse', ...
%     'D:\Full\Processed_PIV_validation\60degree_250light_250hz_1000dt_reverse', ...
    'D:\Full\Processed_PIV_validation\60degree_400light_100hz_3000dt', ...
    'D:\Full\Processed_PIV_validation\60degree_250light_250hz_1000dt', ...
    'D:\Full\Processed_PIV_validation\30degree_400light_100hz_3000dt_reverse', ...
    'D:\Full\Processed_PIV_validation\30degree_250light_250hz_1000dt_reverse', ...
    'D:\Full\Processed_PIV_validation\30degree_400light_100hz_3000dt', ...
    'D:\Full\Processed_PIV_validation\30degree_250light_250hz_1000dt' ...
    };

numBaseDirs = length(base_dirs);

% Define the signals of interest.
    signal_names = { ...
    'Mode 1', ...
    'Mode 2', ...
    'VFlux 1', ...
    'VFlux 11', ...
    'Vortex Center Y1', ...
    'Vortex Center Y2', ...
    'Vortex Size 1', ...
    'Vortex Size 2', ...
    'Convective Velocity 1', ...
    'Convective Velocity 11', ...
    'Delta * 1', ...
    'Delta * 2', ...
    'Rotation Rate 1', ...
    'Rotation Rate 2', ...
    'Momentum Thickness 1', ...
    'Momentum Thickness 2' ...
};


% The cross-correlation computed via xcorr returns a vector of length:
xc_length = 2 * segmentLength - 1;
nan_threshold = 0.3;   % Skip signals with >30% NaN values
%% Determine Sampling Rate
% Extract the sampling rate from the first base directory's name.
[~, baseString, ~] = fileparts(base_dirs{1});
hz_str = regexp(baseString, '\d+hz', 'match');
if ~isempty(hz_str)
    sampling_rate = str2double(hz_str{1}(1:end-2));
else
    warning('Could not determine sampling rate from directory name. Using default 100 Hz.');
    sampling_rate = 100;
end
dt = 1 / sampling_rate;  % Time step in seconds

% Create a lag vector (in samples) for plotting.
lags = -(segmentLength-1):(segmentLength-1);

%% Loop Over Base Directories (Process Each Separately)
for d = 1:numBaseDirs
    base_dir = base_dirs{d};
    fprintf('Processing directory %d of %d:\n   %s\n', d, numBaseDirs, base_dir);
    statistics_dir = fullfile(base_dir, 'Statistics', num2str(N), 'Cam1', 'Instantaneous', 'Calibrated');
    data = load(fullfile(statistics_dir, 'cavity_stats.mat'));
    POD_stats = load(fullfile(statistics_dir, 'Below', 'POD_stats_319x149.mat'));
    time_coefficient = POD_stats.V_svd;

    raw_signals = { ...
        (time_coefficient(1:3600,1)), ... % Mode 1
        (time_coefficient(1:3600,2)), ... % Mode 1
        data.v_flux_all(1:3600, 1), ... % VFlux 1
        data.v_flux_all(1:3600, 11), ... % VFlux 2
        data.vortexCentersY_all(1:3600, 1), ...
        data.vortexCentersY_all(1:3600, 2), ...
        data.vortexSize_all(1:3600, 1), ...
        data.vortexSize_all(1:3600, 2), ...
        data.convective_velocity(1:3600, 1), ...
        data.convective_velocity(1:3600, 11), ...
        data.deltaStar_all(1:3600, 1), ...
        data.deltaStar_all(1:3600, 2), ...
        data.rotation_all(1:3600, 1), ...
        data.rotation_all(1:3600, 2), ...
        data.momentumThickness_all(1:3600, 1), ...
        data.momentumThickness_all(1:3600, 2) ...
    };

    % Check for NaN Threshold and fill missing only if valid
    valid_indices = cellfun(@(sig) mean(isnan(sig)) < nan_threshold, raw_signals);
    signals = cellfun(@(sig) fillmissing(sig, 'linear'), raw_signals(valid_indices), 'UniformOutput', false);
    signals = cell2mat(signals);
    valid_signal_names = signal_names(valid_indices);
    num_signals = width(signals);
    
    % --- Local Accumulators for This Directory ---
    local_xc_sum = zeros(num_signals, num_signals, xc_length);
    local_count_runs = zeros(num_signals, num_signals);
    
    % Construct the statistics directory and load data.

    
    % --- Loop Over Segments ---
    for seg = 1:nSegs
        % Define index range for this segment.
        idx_start = (seg - 1) * segmentLength + 1;
        idx_end   = seg * segmentLength;
        
        raw_signals = { ...
            (time_coefficient(idx_start:idx_end,1)), ... % Mode 1
            (time_coefficient(idx_start:idx_end,2)), ... % Mode 1
            data.v_flux_all(idx_start:idx_end, 1), ... % VFlux 1
            data.v_flux_all(idx_start:idx_end, 11), ... % VFlux 2
            data.vortexCentersY_all(idx_start:idx_end, 1), ...
            data.vortexCentersY_all(idx_start:idx_end, 2), ...
            data.vortexSize_all(idx_start:idx_end, 1), ...
            data.vortexSize_all(idx_start:idx_end, 2), ...
            data.convective_velocity(idx_start:idx_end, 1), ...
            data.convective_velocity(idx_start:idx_end, 11), ...
            data.deltaStar_all(idx_start:idx_end, 1), ...
            data.deltaStar_all(idx_start:idx_end, 2), ...
            data.rotation_all(idx_start:idx_end, 1), ...
            data.rotation_all(idx_start:idx_end, 2), ...
            data.momentumThickness_all(idx_start:idx_end, 1), ...
            data.momentumThickness_all(idx_start:idx_end, 2) ...
        };

        % Check for NaN Threshold and fill missing only if valid
        valid_indices = cellfun(@(sig) mean(isnan(sig)) < nan_threshold, raw_signals);
        signals = cellfun(@(sig) fillmissing(sig, 'linear'), raw_signals(valid_indices), 'UniformOutput', false);
        signals = cell2mat(signals);
        valid_signal_names = signal_names(valid_indices);

        % Skip segment if not enough valid signals
        if size(signals, 2) < 2
            fprintf('Skipping segment %d (too many NaNs)\n', seg);
            continue;
        end
        % Remove the mean from each signal.
        signals = signals - mean(signals, 1);
        
        % Compute cross-correlation for each unique pair (i < j).

        for i = 1:num_signals
            for j = i+1:num_signals
                [xc, ~] = xcorr(signals(:, i), signals(:, j), 'coeff');
                % Accumulate the cross-correlation.
                local_xc_sum(i, j, :) = squeeze(local_xc_sum(i, j, :)) + xc;
                local_count_runs(i, j) = local_count_runs(i, j) + 1;
            end
        end
    end  % End segment loop
    
    % --- Average Over Segments Only ---
    local_avg_xc = zeros(size(local_xc_sum));
    for i = 1:num_signals
        for j = i+1:num_signals
            if local_count_runs(i, j) > 0
                local_avg_xc(i, j, :) = squeeze(local_xc_sum(i, j, :)) / local_count_runs(i, j);
            end
        end
    end
    
    % --- (Optional) Extract Metrics for This Directory ---
    % For example, find the peak correlation, its lag, time shift, and phase difference.
    threshold = 0.15;
    peak_correlations = NaN(num_signals);
    max_lags          = NaN(num_signals);
    time_shifts       = NaN(num_signals);
    phase_differences = NaN(num_signals);
    
    for i = 1:num_signals
        for j = i+1:num_signals
            xc_pair = squeeze(local_avg_xc(i, j, :));
            [~, max_idx] = max(abs(xc_pair));
            peak_val = xc_pair(max_idx);
            
            if abs(peak_val) < threshold
                peak_correlations(i, j) = NaN;
                time_shifts(i, j) = NaN;
                phase_differences(i, j) = NaN;
            else
                peak_correlations(i, j) = peak_val;
                max_lags(i, j) = lags(max_idx);
                time_shifts(i, j) = lags(max_idx) * dt;
                phase_differences(i, j) = time_shifts(i, j) * 360 * sampling_rate;
            end
        end
    end
    
    % --- Plot the Averaged Cross-Correlation Functions ---
    figure('Name', sprintf('Averaged Cross-Correlation Functions for %s', base_dir), 'NumberTitle', 'off');
    plot_idx = 1;
    for i = 1:num_signals
        for j = i+1:num_signals
            subplot(num_signals-1, num_signals-1, plot_idx);
            xc_pair = squeeze(local_avg_xc(i, j, :));
            plot(lags * dt, xc_pair, 'LineWidth', 1.5);
            hold on;
            % Mark the peak correlation.
            [~, max_idx] = max(abs(xc_pair));
            plot(lags(max_idx)*dt, xc_pair(max_idx), 'ro', 'MarkerFaceColor', 'r');
            hold off;
            title(sprintf('%s vs %s', valid_signal_names{i}, valid_signal_names{j}), 'FontSize', 8);
            xlabel('Time Lag (s)', 'FontSize', 7);
            ylabel('Cross-Corr', 'FontSize', 7);
            grid on;
            xlim([-5, 5]);
            ylim([-1, 1]);
            plot_idx = plot_idx + 1;
        end
    end
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    % save figure 
    saveas(gcf, fullfile(statistics_dir, 'Averaged_Cross_Correlation_Functions.fig'));
    saveas(gcf, fullfile(statistics_dir, 'Averaged_Cross_Correlation_Functions.png'));
    
    % Plot Peak Cross-Correlation Heatmap in its own figure
    figure('Name', sprintf('Peak Cross-Correlation for %s', base_dir), 'NumberTitle', 'off');
    h1 = heatmap(valid_signal_names, valid_signal_names, peak_correlations, 'Colormap', parula);
    title('Peak Cross-Correlation');
    h1.ColorLimits = [-1 1];
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    % save figure
    saveas(gcf, fullfile(statistics_dir, 'Peak_Cross_Correlation.fig'));
    saveas(gcf, fullfile(statistics_dir, 'Peak_Cross_Correlation.png'));
    
    % Plot Time Shifts Heatmap in its own figure
    figure('Name', sprintf('Time Shifts for %s', base_dir), 'NumberTitle', 'off');
    h2 = heatmap(valid_signal_names, valid_signal_names, time_shifts, 'Colormap', cool);
    title('Time Shifts (s)');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    % save figure
    saveas(gcf, fullfile(statistics_dir, 'Time_Shifts.fig'));
    saveas(gcf, fullfile(statistics_dir, 'Time_Shifts.png'));


    close all
end  % End directory loop

%% averaged_crosscorr_model.m
% This script processes data from multiple base directories.
% For each directory, the 18,000 data points are split into 5 segments 
% (each of 3,600 points). For each segment the cross‐correlation is computed 
% for a set of signals and then averaged over the segments.
% Plots are produced for each base directory separately.

close all;
clear variables;

%% Parameters and Base Directories
N = 18000;             % Total number of data points per directory
segmentLength = 3600;  % Number of points per segment
nSegs = N / segmentLength;  % Number of segments per directory (should be 5)

base_dirs = { ...
%     'D:\Full\Processed_PIV_validation\90degree_400light_100hz_3000dt', ...
%     'D:\Full\Processed_PIV_validation\90degree_250light_250hz_1000dt', ...
%     'D:\Full\Processed_PIV_validation\60degree_400light_100hz_3000dt_reverse', ...
%     'D:\Full\Processed_PIV_validation\60degree_250light_250hz_1000dt_reverse', ...
    'D:\Full\Processed_PIV_validation\60degree_400light_100hz_3000dt', ...
    'D:\Full\Processed_PIV_validation\60degree_250light_250hz_1000dt', ...
    'D:\Full\Processed_PIV_validation\30degree_400light_100hz_3000dt_reverse', ...
    'D:\Full\Processed_PIV_validation\30degree_250light_250hz_1000dt_reverse', ...
    'D:\Full\Processed_PIV_validation\30degree_400light_100hz_3000dt', ...
    'D:\Full\Processed_PIV_validation\30degree_250light_250hz_1000dt' ...
    };

numBaseDirs = length(base_dirs);

% Define the signals of interest.
    signal_names = { ...
    'Mode 1 Below', ...
    'Mode 2 Below', ...
    'Mode 3 Below', ...
    'Mode 1 Above',...
    'Mode 2 Above',...
    'Mode 3 Above',...
    'Mode 4 Above',...
    'Mode 5 Above',...
    'Mode 6 Above',...
    'Mode 7 Above',...
    'Mode 8 Above',...
    'Mode 9 Above',...
    'Mode 10 Above',...
};


num_signals = length(signal_names);

% The cross-correlation computed via xcorr returns a vector of length:
xc_length = 2 * segmentLength - 1;
nan_threshold = 0.3;   % Skip signals with >30% NaN values
%% Determine Sampling Rate
% Extract the sampling rate from the first base directory's name.
[~, baseString, ~] = fileparts(base_dirs{1});
hz_str = regexp(baseString, '\d+hz', 'match');
if ~isempty(hz_str)
    sampling_rate = str2double(hz_str{1}(1:end-2));
else
    warning('Could not determine sampling rate from directory name. Using default 100 Hz.');
    sampling_rate = 100;
end
dt = 1 / sampling_rate;  % Time step in seconds

% Create a lag vector (in samples) for plotting.
lags = -(segmentLength-1):(segmentLength-1);

%% Loop Over Base Directories (Process Each Separately)
for d = 1:numBaseDirs
    base_dir = base_dirs{d};
    fprintf('Processing directory %d of %d:\n   %s\n', d, numBaseDirs, base_dir);
    statistics_dir = fullfile(base_dir, 'Statistics', num2str(N), 'Cam1', 'Instantaneous', 'Calibrated');
    data = load(fullfile(statistics_dir, 'cavity_stats.mat'));
    POD_stats_below = load(fullfile(statistics_dir, 'Below', 'POD_stats_319x149.mat'));
    POD_stats_above = load(fullfile(statistics_dir, 'Above', 'POD_stats_319x149.mat'));
    

    
    % --- Local Accumulators for This Directory ---
    local_xc_sum = zeros(num_signals, num_signals, xc_length);
    local_count_runs = zeros(num_signals, num_signals);
    
    % Construct the statistics directory and load data.

    
    % --- Loop Over Segments ---
    for seg = 1:nSegs
        % Define index range for this segment.
        idx_start = (seg - 1) * segmentLength + 1;
        idx_end   = seg * segmentLength;
        
        signals = { ...
            POD_stats_below.V_svd(1:3600,1),...
            POD_stats_below.V_svd(1:3600,2),...
            POD_stats_below.V_svd(1:3600,3),...
            POD_stats_above.V_svd(1:3600,1),...
            POD_stats_above.V_svd(1:3600,2),...
            POD_stats_above.V_svd(1:3600,3),...
            POD_stats_above.V_svd(1:3600,4),...
            POD_stats_above.V_svd(1:3600,5),...
            POD_stats_above.V_svd(1:3600,6),...
            POD_stats_above.V_svd(1:3600,7),...
            POD_stats_above.V_svd(1:3600,8),...
            POD_stats_above.V_svd(1:3600,9),...
            POD_stats_above.V_svd(1:3600,10),...  
        };


        valid_signal_names = signal_names;
        signals = cell2mat(signals);

        % Skip segment if not enough valid signals
        if size(signals, 2) < 2
            fprintf('Skipping segment %d (too many NaNs)\n', seg);
            continue;
        end
        % Remove the mean from each signal.
        signals = signals - mean(signals, 1);
        
        % Compute cross-correlation for each unique pair (i < j).

        for i = 1:num_signals
            for j = i+1:num_signals
                [xc, ~] = xcorr(signals(:, i), signals(:, j), 'coeff');
                % Accumulate the cross-correlation.
                local_xc_sum(i, j, :) = squeeze(local_xc_sum(i, j, :)) + xc;
                local_count_runs(i, j) = local_count_runs(i, j) + 1;
            end
        end
    end  % End segment loop
    
    % --- Average Over Segments Only ---
    local_avg_xc = zeros(size(local_xc_sum));
    for i = 1:num_signals
        for j = i+1:num_signals
            if local_count_runs(i, j) > 0
                local_avg_xc(i, j, :) = squeeze(local_xc_sum(i, j, :)) / local_count_runs(i, j);
            end
        end
    end
    
    % --- (Optional) Extract Metrics for This Directory ---
    % For example, find the peak correlation, its lag, time shift, and phase difference.
    threshold = 0.15;
    peak_correlations = NaN(num_signals);
    max_lags          = NaN(num_signals);
    time_shifts       = NaN(num_signals);
    phase_differences = NaN(num_signals);
    
    for i = 1:num_signals
        for j = i+1:num_signals
            xc_pair = squeeze(local_avg_xc(i, j, :));
            [~, max_idx] = max(abs(xc_pair));
            peak_val = xc_pair(max_idx);
            
            if abs(peak_val) < threshold
                peak_correlations(i, j) = NaN;
                time_shifts(i, j) = NaN;
                phase_differences(i, j) = NaN;
            else
                peak_correlations(i, j) = peak_val;
                max_lags(i, j) = lags(max_idx);
                time_shifts(i, j) = lags(max_idx) * dt;
                phase_differences(i, j) = time_shifts(i, j) * 360 * sampling_rate;
            end
        end
    end
    
    % --- Plot the Averaged Cross-Correlation Functions ---
    figure('Name', sprintf('Averaged Cross-Correlation Functions for %s', base_dir), 'NumberTitle', 'off');
    plot_idx = 1;
    for i = 1:num_signals
        for j = i+1:num_signals
            subplot(num_signals-1, num_signals-1, plot_idx);
            xc_pair = squeeze(local_avg_xc(i, j, :));
            plot(lags * dt, xc_pair, 'LineWidth', 1.5);
            hold on;
            % Mark the peak correlation.
            [~, max_idx] = max(abs(xc_pair));
            plot(lags(max_idx)*dt, xc_pair(max_idx), 'ro', 'MarkerFaceColor', 'r');
            hold off;
            title(sprintf('%s vs %s', valid_signal_names{i}, valid_signal_names{j}), 'FontSize', 8);
            xlabel('Time Lag (s)', 'FontSize', 7);
            ylabel('Cross-Corr', 'FontSize', 7);
            grid on;
            xlim([-5, 5]);
            ylim([-1, 1]);
            plot_idx = plot_idx + 1;
        end
    end
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    % save figure 
    saveas(gcf, fullfile(statistics_dir, 'Averaged_MODAL_Cross_Correlation_Functions.fig'));
    saveas(gcf, fullfile(statistics_dir, 'Averaged_MODAL_Cross_Correlation_Functions.png'));
    
    % Plot Peak Cross-Correlation Heatmap in its own figure
    figure('Name', sprintf('Peak Cross-Correlation for %s', base_dir), 'NumberTitle', 'off');
    h1 = heatmap(valid_signal_names, valid_signal_names, peak_correlations, 'Colormap', parula);
    title('Peak Cross-Correlation');
    h1.ColorLimits = [-1 1];
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    % save figure
    saveas(gcf, fullfile(statistics_dir, 'Peak_Cross_MODAL_Correlation.fig'));
    saveas(gcf, fullfile(statistics_dir, 'Peak_Cross_MODAL_Correlation.png'));
    
    % Plot Time Shifts Heatmap in its own figure
    figure('Name', sprintf('Time Shifts for %s', base_dir), 'NumberTitle', 'off');
    h2 = heatmap(valid_signal_names, valid_signal_names, time_shifts, 'Colormap', cool);
    title('Time Shifts (s)');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    % save figure
    saveas(gcf, fullfile(statistics_dir, 'Time_MODAL_Shifts.fig'));
    saveas(gcf, fullfile(statistics_dir, 'Time_MODAL_Shifts.png'));


    close all
end  % End directory loop