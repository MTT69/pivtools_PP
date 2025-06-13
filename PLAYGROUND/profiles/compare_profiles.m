% Clear workspace
clear;
close all;
clc;

% Base directories (same as in Profiles.m)
base_dir = { ...
    'D:\full\Processed_PIV_validation\90degree_400light_100hz_3000dt', ...
    'D:\full\Processed_PIV_validation\90degree_250light_250hz_1000dt', ...
};

% Parameters
N = 18000;
Camera = 1;
window = [16 16];
i = 5; % Window index
x_loc_index = 2; % Use the second x_location (-5) from the original profile analysis

% Create shortened case names
case_names = cell(size(base_dir));
for idx = 1:length(base_dir)
    [~, baseString, ~] = fileparts(base_dir{idx});
    
    % Extract frequency information for the legend
    if contains(baseString, '100hz')
        freq_str = '100 Hz';
    elseif contains(baseString, '250hz')
        freq_str = '250 Hz';
    else
        freq_str = 'Unknown Hz';
    end
    
    case_names{idx} = ['$90^{\circ}$ (' freq_str ')'];
end

% Generate theoretical power law profiles
y_delta_theory = linspace(0, 1.2, 200);
u_uinf_7th = min(y_delta_theory.^(1/7), 1); % 1/7th power law
u_uinf_6th = min(y_delta_theory.^(1/6), 1); % 1/6th power law (changed from 1/5th)
u_uinf_8th = min(y_delta_theory.^(1/8), 1); % 1/8th power law

% Store profiles for comparison plot
all_profiles = cell(length(base_dir), 1);
all_y_values = cell(length(base_dir), 1);
all_deltas = zeros(length(base_dir), 1);

% Process each case individually
for idx = 1:length(base_dir)
    Base = base_dir{idx};
    
    % Load coordinate and velocity data
    calibrated = fullfile(Base, 'CalibratedPIV', num2str(N), ['Cam' num2str(Camera)], 'Instantaneous');
    Co_ords = load(fullfile(calibrated, "Co_ords.mat"));
    
    statistics = fullfile(Base, 'Statistics', num2str(N), ['Cam' num2str(Camera)], 'Instantaneous', 'Calibrated');
    filename = fullfile(statistics, ['MeanStats' num2str(window(1)) 'x' num2str(window(2)) '.mat']);
    meanData = load(filename);
    
    % Extract velocity profile at x = -5
    positive_y_indices = Co_ords.Co_ords(i).y(:,1) > 0;
    x_locations = [-100, -5, 50]; % Same as in original code
    percentile_speed = prctile(meanData.mean_U(positive_y_indices, :), 99, 'all');
    
    % Find x-index for x = -5
    [~, x_idx] = min(abs(Co_ords.Co_ords(i).x(1,:) - x_locations(x_loc_index)));
    
    % Extract profile data
    col_range = max(1, x_idx-2):min(size(Co_ords.Co_ords(i).x, 2), x_idx+2);
    velocity_profile = mean(meanData.mean_U(positive_y_indices, col_range), 2) / percentile_speed;
    y_values = Co_ords.Co_ords(i).y(positive_y_indices, 1);
    
    % Get boundary layer thickness
    [~, delta_idx] = min(abs(velocity_profile - 0.99));
    delta = y_values(delta_idx);
    
    % Store data for comparison plot
    all_profiles{idx} = velocity_profile;
    all_y_values{idx} = y_values;
    all_deltas(idx) = delta;
    
    % Create a new figure for this case
    figure('Position', [100, 100, 900, 700]);
    hold on;
    
    % Plot the actual profile
    plot(y_values/delta, velocity_profile, 'b-o', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', '$\textrm{Measured Profile}$');
    
    % Add theoretical profiles
    plot(y_delta_theory, u_uinf_7th, 'k--', 'LineWidth', 2, 'DisplayName', '$\textrm{1/7th Power Law}$');
    plot(y_delta_theory, u_uinf_6th, 'r:', 'LineWidth', 2, 'DisplayName', '$\textrm{1/6th Power Law}$');
    plot(y_delta_theory, u_uinf_8th, 'g-.', 'LineWidth', 2, 'DisplayName', '$\textrm{1/8th Power Law}$');
    
    % Format plot
    grid on;
    xlabel('$y/\delta$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$u/u_{\infty}$', 'FontSize', 14, 'Interpreter', 'latex');
    title(['$\textrm{Boundary Layer Profile at } x = -5 \textrm{ mm: }$' case_names{idx}], 'FontSize', 16, 'Interpreter', 'latex');
    legend('Location', 'southeast', 'FontSize', 12, 'Interpreter', 'latex');
    axis([0 1.1 0 1.1]);
    set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLabelInterpreter', 'latex');
    
    % Add a text annotation explaining the profiles
    annotation('textbox', [0.15, 0.15, 0.3, 0.12], 'String', ...
        {'$u/u_{\infty} = (y/\delta)^{1/7}$ (dashed line)', ...
         '$u/u_{\infty} = (y/\delta)^{1/6}$ (dotted line)', ...
         '$u/u_{\infty} = (y/\delta)^{1/8}$ (dash-dot line)'}, ...
        'FontSize', 12, 'Interpreter', 'latex', 'BackgroundColor', 'w', 'EdgeColor', 'k');
    
    % Save figure to the case's statistics folder
    saveas(gcf, fullfile(statistics, 'BL-power_law_compare.fig'));
    saveas(gcf, fullfile(statistics, 'BL-power_law_compare.png'));
    
    % Print status
    fprintf('Saved boundary layer comparison for %s\n', case_names{idx});
    
    % Close the figure to avoid memory issues
    close(gcf);
end

% Create comparison subplot figure after processing all cases
if length(base_dir) == 2
    figure('Position', [100, 100, 1200, 600]);
    
    % Add theoretical profiles for reference
    for idx = 1:2
        subplot(1, 2, idx);
        hold on;
        
        % Plot measured profile
        plot(all_y_values{idx}/all_deltas(idx), all_profiles{idx}, 'b-o', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', '$\textrm{Measured Profile}$');
        
        % Plot theoretical profiles
        plot(y_delta_theory, u_uinf_7th, 'k--', 'LineWidth', 1.5, 'DisplayName', '$\textrm{1/7th Power Law}$');
        plot(y_delta_theory, u_uinf_6th, 'r:', 'LineWidth', 1.5, 'DisplayName', '$\textrm{1/6th Power Law}$');
        plot(y_delta_theory, u_uinf_8th, 'g-.', 'LineWidth', 1.5, 'DisplayName', '$\textrm{1/8th Power Law}$');
        
        % Format subplot
        grid on;
        xlabel('$y/\delta$', 'FontSize', 14, 'Interpreter', 'latex');
        ylabel('$u/u_{\infty}$', 'FontSize', 14, 'Interpreter', 'latex');
        title(case_names{idx}, 'FontSize', 16, 'Interpreter', 'latex');
        legend('Location', 'southeast', 'FontSize', 10, 'Interpreter', 'latex');
        axis([0 1.1 0 1.1]);
        set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLabelInterpreter', 'latex');
    end
    
    % Add overall title
    sgtitle('$\textrm{Boundary Layer Profile Comparison at } x = -5 \textrm{ mm}$', 'FontSize', 18, 'Interpreter', 'latex');
    
    % Save comparison figure to working directory
    current_dir = pwd;
    saveas(gcf, fullfile(current_dir, 'BL-profile_comparison.fig'));
    saveas(gcf, fullfile(current_dir, 'BL-profile_comparison.png'));
    
    fprintf('Saved side-by-side profile comparison to %s\n', current_dir);
end

fprintf('Processing complete.\n');