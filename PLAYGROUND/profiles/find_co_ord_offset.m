close all
clear variables

%% Define base directories and parameters
base_dir = { ...
    'D:\full\Processed_PIV_validation\90degree_400light_100hz_3000dt', ...
    'D:\full\Processed_PIV_validation\90degree_250light_250hz_1000dt', ...
    'D:\full\Processed_PIV_validation\60degree_400light_100hz_3000dt_reverse', ...
    'D:\full\Processed_PIV_validation\60degree_250light_250hz_1000dt_reverse', ...
    'D:\full\Processed_PIV_validation\60degree_400light_100hz_3000dt', ...
    'D:\full\Processed_PIV_validation\60degree_250light_250hz_1000dt', ...
    'D:\full\Processed_PIV_validation\30degree_400light_100hz_3000dt_reverse', ...
    'D:\full\Processed_PIV_validation\30degree_250light_250hz_1000dt_reverse', ...
    'D:\full\Processed_PIV_validation\30degree_400light_100hz_3000dt', ...
    'D:\full\Processed_PIV_validation\30degree_250light_250hz_1000dt', ...
};

% Offsets for y_fit adjustment (if needed)
offset = [0, 0.001, 0, 0.001, 0, 0.001, 0, 0.001, 0, 0.001];

scalefactor = 9.53;
N = 18000;
% Settings for data extraction
Camera = 1;
window = [16 16];
i = 5;

for base_idx = 1:length(base_dir)
    Base = base_dir{base_idx};
    [~, baseString, ~] = fileparts(Base);
    
    % Extract dt from the base directory name
    dt_str = regexp(baseString, '\d+dt', 'match');
    dt = str2double(dt_str{1}(1:end-2)) * 1e-6;
    fprintf('Processing file: %s with dt: %.6f\n', baseString, dt);

    statistics = fullfile(Base, 'Statistics', num2str(N), ['Cam' num2str(Camera)], 'Instantaneous', 'Calibrated');
    filename = fullfile(statistics, ['MeanStats' num2str(window(1)) 'x' num2str(window(2)) '.mat']);

end
    