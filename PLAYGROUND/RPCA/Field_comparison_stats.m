% Compare mean estimations
clear variables;
close all;

% Base directory
baseDir = 'D:\full\Processed_PIV_validation\90degree_400light_100hz_3000dt';
statsDir = fullfile(baseDir, 'Statistics\18000\Cam1\Instantaneous\Calibrated\');
calibratedDir = fullfile(baseDir, 'CalibratedPIV\18000\Cam1\');

% Instantaneous files
Inst_Mean = load(fullfile(statsDir, 'MeanStats16x16.mat'));
Inst_Co_ords = load(fullfile(calibratedDir, 'Instantaneous\Co_ords.mat'));

inst_u_prime = Inst_Mean.U_prime_Uprime_mean;
inst_v_prime = Inst_Mean.V_prime_Vprime_mean;


% Ensemble files

Ensemble_Mean = load(fullfile(calibratedDir, 'Ensemble\00001.mat'));
Ensemble_Co_ords = load(fullfile(calibratedDir, 'Ensemble\Co_ords.mat'));


% Add a row to the top of ensemble arrays and remove the bottom row  
%added because the wall location did not line up
ensemble_u = [Ensemble_Mean.piv_result(4).Uturb(1, :); Ensemble_Mean.piv_result(4).Uturb(1:end-1, :)];
ensemble_v = [Ensemble_Mean.piv_result(4).Vturb(1, :); Ensemble_Mean.piv_result(4).Vturb(1:end-1, :)];
% 
% ensemble_u = Ensemble_Mean.piv_result(4).Uturb;
% ensemble_v = Ensemble_Mean.piv_result(4).Vturb;


% Define the colormap for u'
u_values = inst_u_prime(:); % Extract all u' values from instantaneous data
lower_limit_u = prctile(u_values, 1); % 1st percentile for lower limit
upper_limit_u = prctile(u_values, 99); % 99th percentile for upper limit
custommap_u = redbluezero(lower_limit_u, upper_limit_u);

% Define the colormap for v'
v_values = inst_v_prime(:); % Extract all v' values from instantaneous data
lower_limit_v = prctile(v_values, 1); % 1st percentile for lower limit
upper_limit_v = prctile(v_values, 99); % 99th percentile for upper limit
custommap_v = redbluezero(lower_limit_v, upper_limit_v);

% Plot the 2x2 subplot for u' and v' statistics
figure;

% Subplot 1: Instantaneous u'
subplot(2, 2, 1);
imagesc(inst_u_prime);
colormap(custommap_u);
colorbar;
title("Instantaneous u'");
caxis([lower_limit_u, upper_limit_u]); % Use the defined limits for u'
axis image;

% Subplot 2: Instantaneous v'
subplot(2, 2, 2);
imagesc(inst_v_prime);
colormap(custommap_v);
colorbar;
title("Instantaneous v'");
caxis([lower_limit_v, upper_limit_v]); % Use the defined limits for v'
axis image;

% Subplot 3: Ensemble u'
subplot(2, 2, 3);
imagesc(ensemble_u);
colormap(custommap_u); % Use the colormap for u'
colorbar;
title("Ensemble u'");
caxis([lower_limit_u, upper_limit_u]); % Use the same limits as instantaneous u'
axis image;

% Subplot 4: Ensemble v'
subplot(2, 2, 4);
imagesc(ensemble_v);
colormap(custommap_v); % Use the colormap for v'
colorbar;
title("Ensemble v'");
caxis([lower_limit_v, upper_limit_v]); % Use the same limits as instantaneous v'
axis image;


% Calculate multiplicative differences
mult_diff_u = ensemble_u ./ inst_u_prime;
mult_diff_v = ensemble_v ./ inst_v_prime;

% Handle division by zero
mult_diff_u(isinf(mult_diff_u) | isnan(mult_diff_u)) = 0;
mult_diff_v(isinf(mult_diff_v) | isnan(mult_diff_v)) = 0;

% Define color limits (symmetric about 1, as this is a multiplicative factor)
lower_limit_u = -2;
upper_limit_u = 2;

lower_limit_v = -2;
upper_limit_v = 2;

% Create color maps
custommap_u = redbluezero(lower_limit_u, upper_limit_u);
custommap_v = redbluezero(lower_limit_v, upper_limit_v);

% Plot multiplicative differences
figure;

% Multiplicative difference for U'
subplot(2, 1, 1);
imagesc(mult_diff_u, [lower_limit_u, upper_limit_u]);
colormap(custommap_u);
colorbar;
title('Multiplicative Difference in U'' Predictions');
xlabel('X-axis');
ylabel('Y-axis');

% Multiplicative difference for V'
subplot(2, 1, 2);
imagesc(mult_diff_v, [lower_limit_v, upper_limit_v]);
colormap(custommap_v);
colorbar;
title('Multiplicative Difference in V'' Predictions');
xlabel('X-axis');
ylabel('Y-axis');


