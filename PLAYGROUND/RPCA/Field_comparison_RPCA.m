% Compare mean estimations
clear variables;

% Base directory
baseDir = 'D:\full\Processed_PIV\90degree_400light_100hz_3000dt\';
statsDir = fullfile(baseDir, 'Statistics\18000\Cam1\Instantaneous\Calibrated\');
calibratedDir = fullfile(baseDir, 'CalibratedPIV\18000\Cam1\');
uncalibratedDir = fullfile(baseDir, 'UncalibratedPIV\18000\Cam1\');

% Instantaneous files
Inst_Mean = load(fullfile(statsDir, 'MeanStats16x16.mat'));
Inst_Co_ords = load(fullfile(calibratedDir, 'Instantaneous\Co_ords.mat'));
Inst_uncalibrated = load(fullfile(uncalibratedDir, 'Instantaneous\00005.mat'));

% Ensemble files
Ensemble_uncalibrated = load(fullfile(uncalibratedDir, 'Ensemble\00001.mat'));
Ensemble_Mean = load(fullfile(calibratedDir, 'Ensemble\00001.mat'));
Ensemble_Co_ords = load(fullfile(calibratedDir, 'Ensemble\Co_ords.mat'));

% % RPCA files
RPCA_mean = load(fullfile(statsDir, 'RPCA\MeanStats16x16.mat'));

% Extract velocity fields
inst_u = Inst_Mean.mean_U;
inst_v = Inst_Mean.mean_V;

RPCA_u = RPCA_mean.mean_U;
RPCA_v = RPCA_mean.mean_V;

% ensemble_u = Ensemble_Mean.piv_result(4).ux;
% ensemble_v = Ensemble_Mean.piv_result(4).uy;
ensemble_u = [Ensemble_Mean.piv_result(4).ux(1, :); Ensemble_Mean.piv_result(4).ux(1:end-1, :)];
ensemble_v = [Ensemble_Mean.piv_result(4).uy(1, :); Ensemble_Mean.piv_result(4).uy(1:end-1, :)];

% Compute absolute differences for each method
diff_inst_u = abs(ensemble_u - inst_u);
diff_inst_v = abs(ensemble_v - inst_v);

diff_RPCA_u = abs(ensemble_u - RPCA_u);
diff_RPCA_v = abs(ensemble_v - RPCA_v);

% Compute mean absolute error (MAE) for each method
MAE_inst_u = nanmean(diff_inst_u(:));
MAE_inst_v = nanmean(diff_inst_v(:));

MAE_RPCA_u = nanmean(diff_RPCA_u(:));
MAE_RPCA_v = nanmean(diff_RPCA_v(:));

% Define visualization parameters
limit = [-0.1, 0.1];
custommap = redbluezero(limit(1), limit(2));

% Plot difference fields
figure;

subplot(2,2,1);
imagesc(diff_inst_u, limit); colormap(custommap); colorbar;
title('Inst U - Ensemble U Difference');

subplot(2,2,2);
imagesc(diff_inst_v, limit); colormap(custommap); colorbar;
title('Inst V - Ensemble V Difference');

subplot(2,2,3);
imagesc(diff_RPCA_u, limit); colormap(custommap); colorbar;
title('RPCA U - Ensemble U Difference');

subplot(2,2,4);
imagesc(diff_RPCA_v, limit); colormap(custommap); colorbar;
title('RPCA V - Ensemble V Difference');

% Print error metrics to console
fprintf('Mean Absolute Error (U): Inst = %.4f, RPCA = %.4f\n', MAE_inst_u, MAE_RPCA_u);
fprintf('Mean Absolute Error (V): Inst = %.4f, RPCA = %.4f\n', MAE_inst_v, MAE_RPCA_v);

% Define positions for additional profile analysis (example positions)
x_positions = [-105, -20, 0, 20, 105];
y_position = 0; % Modify as needed
