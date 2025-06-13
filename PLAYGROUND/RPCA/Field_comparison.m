% Compare mean estimations
clear variables;

% Base directory
baseDir = 'D:\full\Processed_PIV\90degree_400light_100hz_3000dt';
statsDir = fullfile(baseDir, 'Statistics\18000\Cam1\Instantaneous\Calibrated\');
calibratedDir = fullfile(baseDir, 'CalibratedPIV\18000\Cam1\');
uncalibratedDir = fullfile(baseDir, 'UncalibratedPIV\18000\Cam1\');

% Instantaneous files
Inst_Mean = load(fullfile(statsDir, 'MeanStats16x16.mat'));
Inst_Co_ords = load(fullfile(calibratedDir, 'Instantaneous\Co_ords.mat'));
Inst_uncalibrated = load(fullfile(uncalibratedDir, 'Instantaneous\00005.mat'));

% Ensemble files

Ensemble_Mean = load(fullfile(calibratedDir, 'Ensemble\00001.mat'));
Ensemble_Co_ords = load(fullfile(calibratedDir, 'Ensemble\Co_ords.mat'));

% Extract velocity fields
inst_u = Inst_Mean.mean_U;
inst_v = Inst_Mean.mean_V;

% Add a row to the top of ensemble arrays and remove the bottom row  
%added because the wall location did not line up
ensemble_u = [Ensemble_Mean.piv_result(4).ux(1, :); Ensemble_Mean.piv_result(4).ux(1:end-1, :)];
ensemble_v = [Ensemble_Mean.piv_result(4).uy(1, :); Ensemble_Mean.piv_result(4).uy(1:end-1, :)];
% 
% ensemble_u = Ensemble_Mean.piv_result(4).ux;
% ensemble_v = Ensemble_Mean.piv_result(4).uy;


% Compute absolute differences for each method
diff_inst_u = abs(ensemble_u - inst_u);
diff_inst_v = abs(ensemble_v - inst_v);

% Compute mean absolute error (MAE)
MAE_inst_u = nanmean(diff_inst_u(:));
MAE_inst_v = nanmean(diff_inst_v(:));

% Define visualization parameters
limit = [-0.1, 0.1];
custommap = redbluezero(limit(1), limit(2));

% Plot difference fields
figure;

subplot(1,2,1);
imagesc(diff_inst_u, limit); colormap(custommap); colorbar;
title('Inst U - Ensemble U Difference');
daspect([1 1 1])

subplot(1,2,2);
imagesc(diff_inst_v, limit); colormap(custommap); colorbar;
title('Inst V - Ensemble V Difference');
daspect([1 1 1])

% Print error metrics to console
fprintf('Mean Absolute Error (U): Inst = %.4f\n', MAE_inst_u);
fprintf('Mean Absolute Error (V): Inst = %.4f\n', MAE_inst_v);

% Define positions for additional profile analysis (example positions)
x_positions = [-105, -20, 0, 20, 105];
y_position = 0; % Modify as needed
