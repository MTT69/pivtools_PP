clear; clc;
close all
% Setup paths and parameters
base = 'C:\Users\Lab8-2\Documents\PIVTOOLS\PIVTOOLS\Examples\PIV_data';
run = 7;
synthetic_run = 4;
window = 16;
plot_win = 4;
target_y = 20;  % Desired y position to inspect
type = 'Ensemble';
snr = 10000;

% % Window function parameters
% window_type = 'centered'; % Options: 'hanning', 'hamming', 'tukey', 'gaussian'
% alpha = 0.5;           % Tukey window parameter (0-1), higher = more tapering



% %% 3. Custom Radial Cosine Taper Centered at (9,9)
% % Define coordinate grid
% [x, y] = meshgrid(1:window, 1:window);
% center = [window/2+1, window/2+1];  % Center in MATLAB indexing

% % Compute the Euclidean distance from the center
% r = sqrt((x - center(1)).^2 + (y - center(2)).^2);

% % Set the maximum radius at which the taper goes to zero.
% % Here we choose the maximum distance from (9,9) to the array edge.
% r_max = max([max(center-1), max(window - center)]);

% % Create a radial cosine taper that decays to 0 at r_max
% % The exponent p_radial controls the steepness of the taper.
% p_radial = 3;  % Increase for a faster taper
% window_2d = cos(pi * r / (2 * r_max)).^p_radial;
% % Optionally force taper values to 0 beyond r_max if needed:
% window_2d(r > r_max) = 0;

% % Create a flattened 1D version of the window 
% % using the same reshaping order that MATLAB uses
% window_1d = window_2d(:);


% Load Data
correlation_stats = load(fullfile(base, ['UncalibratedPIV\1000\Cam1\Ensemble\correlation_stats_0000', num2str(run), '.mat']));
Co_ords = load(fullfile('C:\Users\Lab8-2\Documents\PIVTOOLS\PIVTOOLS\Examples\Synthetic_profiles\Co_ords.mat'));

% Optional: Mask out the top half if needed
[m, n] = size(Co_ords.Co_ords(synthetic_run).x);
midpoint = floor(m / 2);
Co_ords.Co_ords(synthetic_run).y_plus(1:midpoint, :) = NaN;

% Find the closest y location
y_plus_flat = Co_ords.Co_ords(synthetic_run).y_plus(:);
[~, closest_idx] = min(abs(y_plus_flat - target_y));

% Get the (row, col) index of this closest y
[row_y, col_y] = ind2sub(size(Co_ords.Co_ords(synthetic_run).y_plus), closest_idx);
actual_y = Co_ords.Co_ords(synthetic_run).y_plus(row_y, col_y);
fprintf('Target y+ = %.2f, Closest y+ found = %.2f at (row = %d)\window', target_y, actual_y, row_y);

% Prepare available x positions along this y row
available_cols = 1:window;  % All columns in this row are potential x positions

% Randomly select plot_count x positions
plot_count = 1;
selected_cols = available_cols(randperm(length(available_cols), plot_count));

% Create a figure with subplots
for i = 1:plot_count
    idx_x = selected_cols(i);   % Current x index
    idx_y = row_y;              % Fixed y index (closest y row)

    % Extract original correlation vectors
    corr_flat = correlation_stats.correlation_plane(:, idx_y, idx_x);
    AA_flat = correlation_stats.pointspreadA(:, idx_y, idx_x);
    BB_flat = correlation_stats.pointspreadB(:, idx_y, idx_x);
    
    % Find the peak height in each correlation map
    peak_corr = max(corr_flat);
    peak_AA = max(AA_flat);
    peak_BB = max(BB_flat);
    
    % Apply 3% threshold filter (set values below 3% of peak to zero)
    corr_flat_windowed = corr_flat;
    AA_flat_windowed = AA_flat;
    BB_flat_windowed = BB_flat;
    
    corr_flat_windowed(corr_flat < 0.03 * peak_corr) = 0;
    AA_flat_windowed(AA_flat < 0.03 * peak_AA) = 0;
    BB_flat_windowed(BB_flat < 0.03 * peak_BB) = 0;
    

    % Extract original and windowed correlation maps
    corr_map = reshape(corr_flat, window, window);
    AA_map = reshape(AA_flat, window, window);
    BB_map = reshape(BB_flat, window, window);
    
    % Reshape windowed 1D vectors to 2D maps
    corr_map_windowed = reshape(corr_flat_windowed, window, window);
    AA_map_windowed = reshape(AA_flat_windowed, window, window);
    BB_map_windowed = reshape(BB_flat_windowed, window, window);
    
    % Row 1: Original maps
    % Auto-correlation A
    subplot(plot_count*2, 3, (i-1)*6 + 1);
    surf(AA_map);
    title(sprintf('Original Auto-A (x=%d, y=%d)', idx_x, idx_y));
    shading interp;
    %view(45,35);
    colorbar;
    
    % Auto-correlation B
    subplot(plot_count*2, 3, (i-1)*6 + 2);
    surf(BB_map);
    title(sprintf('Original Auto-B (x=%d, y=%d)', idx_x, idx_y));
    shading interp;
    %view(45,35);
    colorbar;
    
    % Cross-correlation
    subplot(plot_count*2, 3, (i-1)*6 + 3);
    surf(corr_map);
    title(sprintf('Original Cross-correlation (x=%d, y=%d)', idx_x, idx_y));
    shading interp;
    %view(45,35);
    colorbar;
    
    % Row 2: Windowed maps
    % Auto-correlation A
    subplot(plot_count*2, 3, (i-1)*6 + 4);
    surf(AA_map_windowed);
    title(sprintf('Windowed '));
    shading interp;
    %view(45,35);
    colorbar;
    
    % Auto-correlation B
    subplot(plot_count*2, 3, (i-1)*6 + 5);
    surf(BB_map_windowed);
    title(sprintf('Windowed'));
    shading interp;
    %view(45,35);
    colorbar;
    
    % Cross-correlation
    subplot(plot_count*2, 3, (i-1)*6 + 6);
    surf(corr_map_windowed);
    title(sprintf('Windowed'));
    shading interp;
    %view(45,35);
    colorbar;
end

% Add a super title
sgtitle(sprintf('Correlation Maps at y+ ≈ %.2f with Windowing', actual_y), 'FontSize', 14, 'FontWeight', 'bold');
