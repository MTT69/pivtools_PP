close all
clear variables
base_dir = { ...
    'D:\full\Processed_PIV_validation\90degree_400light_100hz_3000dt', ...
%     'D:\full\Processed_PIV_validation\90degree_250light_250hz_1000dt', ...
%     'D:\full\Processed_PIV_validation\60degree_400light_100hz_3000dt_reverse', ...
%     'D:\full\Processed_PIV_validation\60degree_250light_250hz_1000dt_reverse', ...
%     'D:\full\Processed_PIV_validation\60degree_400light_100hz_3000dt', ...
%     'D:\full\Processed_PIV_validation\60degree_250light_250hz_1000dt', ...
%     'D:\full\Processed_PIV_validation\30degree_400light_100hz_3000dt_reverse', ...
%     'D:\full\Processed_PIV_validation\30degree_250light_250hz_1000dt_reverse', ...
%     'D:\full\Processed_PIV_validation\30degree_400light_100hz_3000dt', ...
%     'D:\full\Processed_PIV_validation\30degree_250light_250hz_1000dt', ...
};

scalefactor = 9.53;
N=18000;

for base_idx = 1:length(base_dir)
    disp(['Processing base directory ', num2str(base_idx), ' of ', num2str(length(base_dir))]);
    Base = base_dir{base_idx};
    % Extract the base string after the last backslash
    [~, baseString, ext] = fileparts(Base);
    Camera = 1;
    window = [16 16];
    i = 5;
    lower_limit = -0.25;
    upper_limit = 0.25;
    statistics = fullfile(Base, 'Statistics', num2str(N), ['Cam' num2str(Camera)], 'Instantaneous', 'Calibrated','Below');
    load(fullfile(statistics,'POD_stats_319x149.mat'))
    
    % Extract the sampling rate from the base directory name
    hz_str = regexp(baseString, '\d+hz', 'match');
    sampling_rate = str2double(hz_str{1}(1:end-2));
    
    % Initialize structure to store results
    results = struct();
    
    for mode = 1:3
        coeff = S_svd(mode,mode)* (V_svd(:,mode));

        xc_length = 2 * 3600 - 1;
    
        % Initialize an accumulator for the autocorrelation and a counter.
        local_ac_sum = zeros(1, xc_length);
        
        for seg = 1:5
            % Define index range for this segment.
            idx_start = (seg - 1) * 3600 + 1;
            idx_end   = seg * 3600;
            
            % Extract the segment from the data.
            seg_data = coeff(idx_start:idx_end)';
            
            
            % Compute autocorrelation for this segment.
            ac_seg = xcorr(seg_data, 'unbiased');
            
            % Accumulate the autocorrelation.
            local_ac_sum = local_ac_sum + ac_seg;
        end
        autocorr = local_ac_sum / 5;
        
        % Plot the autocorrelation
        figure;
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
        plot(autocorr);
        title(['Unbiased Autocorrelation for Mode ', num2str(mode)]);
        xlabel('Lag');
        ylabel('Autocorrelation');
        xlim([3100,4100])
        
        % Manually select the primary and secondary peaks
        [x, y] = ginput(2);
        locs = round(x);
        pks = y;
        
        % Check if the user pressed escape (selected points are too close)
        if length(locs) < 2 || abs(locs(2) - locs(1)) < 10
            locs = nan;
            pks = nan;
        end
        
        if length(locs) > 1
            period_lag = locs(2) - locs(1);
            period_time = period_lag / sampling_rate;
            frequency = 1 / period_time;
        else
            period_lag = NaN; % Not enough peaks to determine the period
            period_time = NaN;
            frequency = NaN;
        end
        
        % Store results in the structure
        results(mode).autocorr = autocorr;
        results(mode).peak_locations = locs;
        results(mode).period_lag = period_lag;
        results(mode).period_time = period_time;
        results(mode).frequency = frequency;
        
        % Update the plot with selected peaks and save the figure
        hold on;
        plot(locs, pks, 'ro'); % Mark the peaks
        hold off;
        title(['Unbiased Autocorrelation for Mode ', num2str(mode), ...
               ' - PerioF: ', num2str(period_time, '%.2f'), ' s', ...
               ' - Frequency: ', num2str(frequency, '%.2f'), ' Hz']);
        
        
        saveas(gcf, fullfile(statistics, ['autocorr_mode_', num2str(mode), '.png']));
        saveas(gcf, fullfile(statistics, ['autocorr_mode_', num2str(mode), '.fig']));
        close gcf
    end
    close all
    % Save the results structure
    save(fullfile(statistics, 'autocorr_results.mat'), 'results');

    num_modes = 3; % First three modes
    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.5, 1, 0.4]); % Adjust figure size for horizontal layout
    
    for mode = 1:num_modes
        % Extract autocorrelation for mode
        autocorr_values = results(mode).autocorr;
        
        % Create subplot
        subplot(1, num_modes, mode);
        plot(autocorr_values, 'b'); % Plot the autocorrelation function
        xlabel('Lag');
        ylabel('Autocorrelation');
        
        % Add title based on previously generated info
        title(['Mode ', num2str(mode), ' - Freq: ', num2str(results(mode).frequency, '%.2f'), ' Hz']);
        
        % Adjust formatting
        xlim([3100,4100])
        grid on;
    end
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    % Save the figure
    saveas(gcf, fullfile(statistics, 'autocorrelation_modes_1_3.jpg'));
    disp('Autocorrelation subplot saved as JPEG.');
    close(gcf)


    
    disp('stop')
    
end

