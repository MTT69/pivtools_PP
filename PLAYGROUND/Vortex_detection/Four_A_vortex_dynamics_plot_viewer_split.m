close all;
clear variables;

% (Set N to the same value you used when processing the results.)
N = 18000;
CameraNo =1;

% List of base directories (adjust the paths as needed)
base_dirs = { ...
    'D:\full\Processed_PIV_validation\90degree_400light_100hz_3000dt', ...
    'D:\full\Processed_PIV_validation\90degree_250light_250hz_1000dt', ...
    'D:\full\Processed_PIV_validation\60degree_400light_100hz_3000dt_reverse', ...
    'D:\full\Processed_PIV_validation\60degree_250light_250hz_1000dt_reverse', ...
    'D:\full\Processed_PIV_validation\60degree_400light_100hz_3000dt', ...
    'D:\full\Processed_PIV_validation\60degree_250light_250hz_1000dt', ...
    'D:\full\Processed_PIV_validation\30degree_400light_100hz_3000dt_reverse', ...
    'D:\full\Processed_PIV_validation\30degree_250light_250hz_1000dt_reverse', ...
};


% List of figure filenames to process
fig_files = { ...
    'Autocorrelation_vflux1.fig', ...
    'Autocorrelation_vflux11.fig', ...
    'Autocorrelation_vortexCenterY_1.fig', ...
    'Autocorrelation_vortexCenterY_2.fig', ...
    'Autocorrelation_vortexSize1.fig',...
    'Autocorrelation_vortexSize2.fig',...
    'Autocorrelation_convective velocity1.fig',...
    'Autocorrelation_convective velocity11.fig',...
    'Autocorrelation_delta_star1.fig',...
    'Autocorrelation_delta_star2.fig',...
    'Autocorrelation_rotationRate1.fig',...
    'Autocorrelation_rotationRate2.fig',...
    'Autocorrelation_Momentum thickness1.fig',...
    'Autocorrelation_Momentum thickness1.fig',...
};

% Loop over each base directory
for base_idx = 1:length(base_dirs)
    base_dir = base_dirs{base_idx};
    disp(['Processing base directory ' num2str(base_idx) ' of ' num2str(length(base_dirs))]);
    statistics_dir = fullfile(base_dir, 'Statistics', num2str(N), 'Cam1', 'Instantaneous', 'Calibrated');
    data = load(fullfile(statistics_dir,'cavity_stats.mat'));
    disp(['Vortex 1 rotation rate rev/s: ', num2str(mean(data.rotation_all(:,1), 'omitnan'))])
    disp(['Vortex 2 rotation rate rev/s: ', num2str(mean(data.rotation_all(:,2), 'omitnan'))])
    
    [~, baseString, ~] = fileparts(base_dir);
    hz_str = regexp(baseString, '\d+hz', 'match');
    if ~isempty(hz_str)
        sampling_rate = str2double(hz_str{1}(1:end-2));
    else
        warning('Could not determine sampling rate from directory name. Using default sampling rate 100 Hz.');
        sampling_rate = 100;
    end
    
    frequencies = struct();
    frequencies.names = fig_files;
    
    for fig_idx = 1:length(fig_files)
        fig_path = fullfile(statistics_dir, fig_files{fig_idx});
        if exist(fig_path, 'file')
            fig_handle = openfig(fig_path, 'new');
            ax_handles = findall(fig_handle, 'Type', 'axes');
            for a = 1:length(ax_handles)
                ax = ax_handles(a);
                xlim(ax, [3100, 4100]); % Set x-axis limits
                titleHandle = get(ax, 'Title');
                titleStr = get(titleHandle, 'String');
                
                if contains(lower(titleStr), 'autocorrelation')
                    line_handle = findobj(ax, 'Type', 'line');
                    if ~isempty(line_handle)
                        lh = line_handle(1);
                        xData = get(lh, 'XData');
                        yData = get(lh, 'YData');
                        
                        idx_range = (xData >= 3100) & (xData <= 3600);
                        if any(idx_range)
                            xData_range = xData(idx_range);
                            yData_range = yData(idx_range);
                            
                            [x, y] = ginput(2);
                            locs = round(x);
                            pks = y;
                            
                            if length(locs) < 2 || abs(locs(2) - locs(1)) < 10
                                locs = nan;
                                pks = nan;
                            end
                            
                            if length(locs) >= 2
                                period_lag = locs(2) - locs(1);
                                period_time = period_lag / sampling_rate;
                                frequency = 1 / period_time;
                                disp(frequency)
                            else
                                period_time = NaN;
                                frequency = NaN;
                            end
                            frequencies(fig_idx).freq= frequency;
                            newTitle = sprintf('%s - Period: %.2f s - Frequency: %.2f Hz', fig_files{fig_idx}, period_time, frequency);
                            set(titleHandle, 'String', newTitle, 'FontSize', 30);
                        else
                            set(titleHandle, 'FontSize', 30);
                        end
                    else
                        set(titleHandle, 'FontSize', 30);
                    end
                else
                    set(titleHandle, 'FontSize', 16);
                end
            end
            
            saveas(fig_handle, fig_path);
            save_dir = fullfile(base_dirs{base_idx}, 'Statistics', num2str(N), ['Cam', num2str(CameraNo)], 'Instantaneous', 'Calibrated');
            save_path = fullfile(save_dir, 'Vortex_internal_frequencies.mat');
            save(save_path,"frequencies")
            disp(['Updated and saved figure: ' fig_files{fig_idx}]);
            disp(['frequency' num2str(frequency) 'hz'])
            figure(1);
            close(fig_handle);
        else
            warning(['File not found: ' fig_path]);
        end
    end
end
