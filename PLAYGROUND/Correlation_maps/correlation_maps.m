clear; clc;
close all
% Setup paths and parameters
base = 'C:\Users\Lab8-2\Documents\PIVTOOLS\PIVTOOLS\Examples\PIV_data';
run = 8;
synthetic_run = 5;
window = 32;
plot_win = 32;
target_y = 0;  % Desired y position to inspect
type = 'Ensemble';
snr = 5;

% Load Data
correlation_stats = load(fullfile(base, ['UncalibratedPIV\1000\Cam1\Ensemble\correlation_stats_0000', num2str(run), '.mat']));
piv_result = load(fullfile(base, 'UncalibratedPIV\1000\Cam1\Ensemble\00001.mat'));
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
fprintf('Target y+ = %.2f, Closest y+ found = %.2f at (row = %d)\n', target_y, actual_y, row_y);

% Prepare available x positions along this y row
available_cols = 1:n;  % All columns in this row are potential x positions

% Randomly select 5 x positions
selected_cols = available_cols(randperm(length(available_cols), 5));

% Create a figure with subplots

for i = 1:5
    idx_x = selected_cols(i);   % Current x index
    idx_y = row_y;              % Fixed y index (closest y row)

    % Extract correlation map - adjust this if your data indexing differs
    % Assumes correlation_plane(:, x_idx, y_idx) shape
    corr_map = reshape(correlation_stats.correlation_plane(:, idx_y, idx_x), window, window);
    AA_map = reshape(correlation_stats.pointspreadA(:, idx_y, idx_x), window, window);
    BB_map = reshape(correlation_stats.pointspreadB(:, idx_y, idx_x), window, window);
    
    % Extract ellipse parameters
    sig_x = piv_result.piv_result(run).sig_AB_x(idx_y, idx_x);
    sig_y = piv_result.piv_result(run).sig_AB_y(idx_y, idx_x);
    rotation = piv_result.piv_result(run).rotation_PD(idx_y, idx_x);
    
    % Extract predictor fields and displacement vectors
    predx = piv_result.piv_result(run).Predictor_Field(idx_y, idx_x, 2);
    predy = piv_result.piv_result(run).Predictor_Field(idx_y, idx_x, 1); % Corrected from Predcitor_Field
    ux = piv_result.piv_result(run).ux(idx_y, idx_x);
    uy = piv_result.piv_result(run).uy(idx_y, idx_x);
    
%     % Extract AA parameters (also used for BB)
%     sig_A_x = piv_result.piv_result(run).sig_A_x(idx_y, idx_x);
%     sig_A_y = piv_result.piv_result(run).sig_A_y(idx_y, idx_x);
%     sig_A_xy = piv_result.piv_result(run).sig_A_xy(idx_y, idx_x);
%     
%     % Calculate rotation using eigendecomposition
%     rotation_A = calculate_rotation(sig_A_x, sig_A_y, sig_A_xy);
%     
%     % Get eigenvalues (calculated from covariance matrix)
%     [eigenA_x, eigenA_y] = calculate_eigenvalues(sig_A_x, sig_A_y, sig_A_xy);
% 
%     % Plot the correlation map
%     subplot(3, 5, i);
%     imagesc(corr_map);
%     colorbar;
%     axis image;
%     title(sprintf('Corr Map %d (x idx = %d)', i, idx_x));
    hold on;
    
%     % Draw ellipse on correlation map
%     if ~isnan(sig_x) && ~isnan(sig_y) && ~isnan(rotation)
%         % Define ellipse center with corrected position for correlation map
%         center_x = window/2 + 1 + (ux - predx);
%         center_y = window/2 + 1 + (uy - predy);
%         
%         % Scale factors to make ellipse visible (adjust as needed)
%         scale = 1.5;
%         
%         % Draw rotated ellipse
%         draw_ellipse(center_x, center_y, scale*sig_x, scale*sig_y, rotation);
%         
%         % Add text annotation with values
%         text(1, 1, sprintf('x: %.2f, y: %.2f\nrot: %.1f°', sig_x, sig_y, rotation), ...
%             'Color', 'white', 'FontSize', 8, 'VerticalAlignment', 'top', 'BackgroundColor', [0 0 0 0.5]);
%     end
    
    % Extract and display Uturb and Vturb values
    if isfield(piv_result.piv_result(run), 'Uturb') && isfield(piv_result.piv_result(run), 'Vturb')
        Uturb_val = piv_result.piv_result(run).Uturb(idx_y, idx_x);
        Vturb_val = piv_result.piv_result(run).Vturb(idx_y, idx_x);
        text(0.5, -0.3, sprintf('Uturb: %.3f, Vturb: %.3f', Uturb_val, Vturb_val), ...
            'Color', 'black', 'FontSize', 9, 'HorizontalAlignment', 'center', 'Units', 'normalized');
    end
    hold off;

    % Plot the AA map
    subplot(3, 5, i + 5);
    imagesc(AA_map);
    colorbar;
    axis image;
    title(sprintf('AA Map %d (x idx = %d)', i, idx_x));
    hold on;
    
    % Draw ellipse on AA map
    if ~isnan(sig_A_x) && ~isnan(sig_A_y) && ~isnan(sig_A_xy)
        center_x = window/2 + 1;
        center_y = window/2 + 1;
        scale = 1.5;
        
        % Draw rotated ellipse using eigenvalues
        draw_ellipse(center_x, center_y, scale*sqrt(eigenA_x), scale*sqrt(eigenA_y), rotation_A);
        
        % Add text annotation with values
        text(1, 1, sprintf('x: %.2f, y: %.2f\nrot: %.1f°', sqrt(eigenA_x), sqrt(eigenA_y), rotation_A), ...
            'Color', 'white', 'FontSize', 8, 'VerticalAlignment', 'top', 'BackgroundColor', [0 0 0 0.5]);
    end
    hold off;

    % Plot the BB map
    subplot(3, 5, i + 10);
    imagesc(BB_map);
    colorbar;
    axis image;
    title(sprintf('BB Map %d (x idx = %d)', i, idx_x));
    hold on;
    
    % Draw ellipse on BB map - using same sig_A parameters
    if ~isnan(sig_A_x) && ~isnan(sig_A_y) && ~isnan(sig_A_xy)
        center_x = window/2 + 1;
        center_y = window/2 + 1;
        scale = 1.5;
        
        % Draw rotated ellipse using same eigenvalues as AA
        draw_ellipse(center_x, center_y, scale*sqrt(eigenA_x), scale*sqrt(eigenA_y), rotation_A);
        
        % Add text annotation with values
        text(1, 1, sprintf('x: %.2f, y: %.2f\nrot: %.1f°', sqrt(eigenA_x), sqrt(eigenA_y), rotation_A), ...
            'Color', 'white', 'FontSize', 8, 'VerticalAlignment', 'top', 'BackgroundColor', [0 0 0 0.5]);
    end
    hold off;
end

% Overall title
sgtitle(sprintf('y+ = %.2f (Actual y = %.2f), Window = %d, Type = %s SNR = %s', target_y, actual_y, plot_win, type, num2str(snr)));

% Helper function to draw rotated ellipse
function h = draw_ellipse(x, y, a, b, angle_degrees)
    % Convert angle to radians
    angle_rad = deg2rad(angle_degrees);
    
    % Generate points for an ellipse
    t = linspace(0, 2*pi, 100);
    X = a * cos(t);
    Y = b * sin(t);
    
    % Rotation matrix
    R = [cos(angle_rad), -sin(angle_rad); sin(angle_rad), cos(angle_rad)];
    
    % Rotate points
    points = R * [X; Y];
    
    % Translate to center
    X_rot = points(1, :) + x;
    Y_rot = points(2, :) + y;
    
    % Plot ellipse
    h = plot(X_rot, Y_rot, 'r-', 'LineWidth', 2);
end

% Function to calculate rotation angle using eigendecomposition
function rotation = calculate_rotation(sx, sy, sxy)
    % Create covariance matrix
    covar = [sx, sxy; sxy, sy];
    
    % Compute eigenvalues and eigenvectors
    [eigenvectors, eigenvalues] = eig(covar);
    eigenvalues = diag(eigenvalues);
    
    % Calculate the angles of the eigenvectors
    angles = atan2(eigenvectors(2, :), eigenvectors(1, :));
    angles = rad2deg(angles);
    
    % Determine the correct angle based on the conditions
    if eigenvectors(1, 1) < 0
        angle1 = 180 - abs(angles(1));
    else
        angle1 = abs(angles(1));
    end
    
    if eigenvectors(2, 1) < 0
        angle2 = 180 - abs(angles(2));
    else
        angle2 = abs(angles(2));
    end
    
    % Assign correct angle
    if angle1 < angle2
        rotation = angle1;
    else
        rotation = angle2;
    end
end

% Function to calculate eigenvalues from covariance parameters
function [eigenx, eigeny] = calculate_eigenvalues(sx, sy, sxy)
    % Create covariance matrix
    covar = [sx, sxy; sxy, sy];
    
    % Compute eigenvalues
    eigenvalues = eig(covar);
    
    % Sort eigenvalues (largest first)
    eigenvalues = sort(eigenvalues, 'descend');
    
    % Return the eigenvalues
    eigenx = eigenvalues(1);
    eigeny = eigenvalues(2);
end

