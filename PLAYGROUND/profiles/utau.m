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
offset = [0, 0.00, 0, 0.001, 0, 0.001, 0, 0.001, 0, 0.001];

scalefactor = 9.53;
N = 18000;

% There are 10 directories, which we assume come in pairs (5 pairs total)
nPairs = length(base_dir) / 2;

% Preallocate a structure array to store processing results for each base directory.
results = struct();

% We'll use these colors to distinguish the two members of each pair.
colors = {'b','r'};

%% Process each base directory and store results
for base_idx = 1:length(base_dir)
    Base = base_dir{base_idx};
    [~, baseString, ~] = fileparts(Base);
    
    % Extract dt from the base directory name
    dt_str = regexp(baseString, '\d+dt', 'match');
    dt = str2double(dt_str{1}(1:end-2)) * 1e-6;
    fprintf('Processing file: %s with dt: %.6f\n', baseString, dt);
    
    % Settings for data extraction
    Camera = 1;
    window = [16 16];
    i = 5;
    
    % Load the necessary data files
    calibrated = fullfile(Base, 'CalibratedPIV', num2str(N), ['Cam' num2str(Camera)], 'Instantaneous');
    Co_ords = load(fullfile(calibrated, "Co_ords.mat"));
    VelData = load(fullfile(calibrated, "00001.mat"));
    
    statistics = fullfile(Base, 'Statistics', num2str(N), ['Cam' num2str(Camera)], 'Instantaneous', 'Calibrated');
    filename = fullfile(statistics, ['MeanStats' num2str(window(1)) 'x' num2str(window(2)) '.mat']);
    meanData = load(filename);
    
    % Use only the y > 0 region
    positive_y_indices = Co_ords.Co_ords(i).y(:,1) > 0;
    y = Co_ords.Co_ords(i).y(positive_y_indices, 1);
    
    % Determine U_inf from the 99th percentile speed (all columns)
    percentile_speed = prctile(meanData.mean_U(positive_y_indices, :), 99, 'all');
    disp(["U-inf", percentile_speed])
    
    % Select x-location(s) for profile extraction
    x_locations = [-5];
    x_index = cell(1, numel(x_locations));
    velocity_profiles = cell(numel(x_locations),1);
    for idx = 1:length(x_locations)
        x_loc = x_locations(idx);
        [~, x_index{idx}] = min(abs(Co_ords.Co_ords(i).x(1,:) - x_loc));
        % Average over a few columns around the chosen x-location
        col_range = max(1, x_index{idx}-2):min(size(Co_ords.Co_ords(i).x,2), x_index{idx}+2);
        velocity_profiles{idx} = mean(meanData.mean_U(positive_y_indices, col_range), 2);
    end
    
    % Compute momentum thickness (theta)
    momentum_thickness = zeros(length(x_locations), 1);
    for idx = 1:length(x_locations)
        U = velocity_profiles{idx}/percentile_speed; % velocity profile at current x-location
        % Here we assume U_inf = 1 for the integration; we then normalize using percentile_speed
        U_inf = 1; 
        % Find the index where U is closest to U_inf (i.e. freestream)
        [~, y_delta_index] = min(abs(U - U_inf));
        y_delta = y(y_delta_index); % y-coordinate corresponding to U_inf
        
        % Find the index of the y-coordinate closest to zero
        [~, y_zero_index] = min(abs(y));
        
        % Restrict y and U to the range [y_zero, y_delta]
        y_restricted = flipud(y(y_delta_index: y_zero_index));
        U_restricted = flipud(U(y_delta_index: y_zero_index));
        
        % Calculate the integrand for momentum thickness
        integrand = (U_restricted ./ U_inf) .* (1 - (U_restricted ./ U_inf));
        theta = trapz(y_restricted, integrand);
        momentum_thickness(idx) = theta;
    end
    % Calculate Re_theta (using the first x_location; note conversion of theta from mm to m)
    nu = 1e-6; % kinematic viscosity
    Re_theta = percentile_speed * momentum_thickness(1)*1e-3 / nu;
    
    fprintf('Re_theta @ x = %d --> %.4e\n', x_locations(1), Re_theta);
    
    % Determine delta from the y-location where U reaches the 99th percentile
    [~, y_percentile_index] = min(abs(meanData.mean_U(positive_y_indices, :) - percentile_speed), [], 1);
    y_percentile_position = Co_ords.Co_ords(i).y(positive_y_indices, 1);
    y_percentile_position = y_percentile_position(y_percentile_index);
    
    fprintf('delta: %.4f mm\n', y_percentile_position(x_index{1}));
    fprintf('U_inf: %.4f m/s\n', percentile_speed);
    
    %% --- Inner Scaling (Log-Law) Processing ---
    % Fixed parameters for the log-law model
    kappa = 0.37;
    B = 4.3;
    nu = 1e-6;
    
    % Define fitting region (from 5% to 15% of delta)
    y_min = 0.05 * y_percentile_position(x_index{1});
    y_max = 0.15 * y_percentile_position(x_index{1});
    fit_indices = (y >= y_min) & (y <= y_max);
    
    % Convert y to meters and apply offset if needed
    y_fit = y(fit_indices) * 1e-3;
    y_fit = y_fit - offset(base_idx);
    U_fit = meanData.mean_U(fit_indices, x_index{1});
    
    % Initial guess for u_tau
    u_tau_initial_guess = 0.05 * percentile_speed;
    
    % Define the log-law function (in physical units)
    log_law = @(u_tau, y_fit) ((u_tau / kappa) .* log(y_fit .* u_tau / nu) + B * u_tau);
    
    % Optimization options
    lsqoptions = optimoptions('lsqcurvefit', ...
        'Display', 'off', ...
        'Algorithm', 'levenberg-marquardt', ...
        'TolFun', 1e-10, ...
        'TolX', 1e-10, ...
        'MaxIter', 1000, ...
        'MaxFunEvals', 2000);
    
    % Fit for u_tau using lsqcurvefit
    [u_tau, ~, ~, ~, ~] = lsqcurvefit(log_law, u_tau_initial_guess, y_fit, U_fit, [], [], lsqoptions);
    
    % Convert experimental data to inner units
    y_plus_exp = (y_fit * u_tau) / nu;
    U_plus_exp = U_fit / u_tau;
    
    % Compute the theoretical inner scaling curve for a range of y+ values
    y_plus_th = logspace(0.1, 3, 100);
    U_plus_th = (1/kappa) * log(y_plus_th) + B;
    
    % Store inner scaling results in the results structure
    results(base_idx).baseString   = baseString;
    results(base_idx).y_plus_th    = y_plus_th;
    results(base_idx).U_plus_th    = U_plus_th;
    results(base_idx).y_plus_exp   = y_plus_exp;
    results(base_idx).U_plus_exp   = U_plus_exp;
    results(base_idx).u_tau        = u_tau;
    
    %% --- Outer Scaling Processing ---
    % Outer scaling uses the full profile at the chosen x-location.
    y_over_delta = y(positive_y_indices) ./ y_percentile_position(x_index{1});
    U_over_Uinf  = velocity_profiles{1} ./ percentile_speed;
    
    results(base_idx).y_over_delta = y_over_delta;
    results(base_idx).U_over_Uinf  = U_over_Uinf;
    results(base_idx).percentile_speed = percentile_speed;
    
    % Also store additional stats for annotation:
    results(base_idx).theta = momentum_thickness(1)*1e-3;  % convert mm to m
    results(base_idx).Re_theta = Re_theta;
    results(base_idx).delta = y_percentile_position(x_index{1}); % in mm
    results(base_idx).Uinf = percentile_speed;
end

%% --- Plotting: Create 5 separate figures (each with a 1x2 subplot layout) ---
for pair = 1:nPairs
    % Create a new figure for this pair
    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    
    %% Left subplot: Inner Scaling (Log-Law)
    subplot(1,2,1);
    hold on;
    for j = 1:2
        idx = 2*(pair-1) + j;
        % Plot the theoretical log-law (dashed) for this case
        semilogx(results(idx).y_plus_th, results(idx).U_plus_th, [colors{j} '--'], 'LineWidth', 1.5, 'HandleVisibility','off');
        % Plot the experimental inner-scaled data (markers)
        scatter(results(idx).y_plus_exp, results(idx).U_plus_exp, 60, colors{j}, 'filled', ...
            'DisplayName', results(idx).baseString);
    end
    grid on;
    xlabel('$y^+$','Interpreter','latex','FontSize',12);
    ylabel('$U^+$','Interpreter','latex','FontSize',12);
    title(sprintf('Inner Scaling (Pair %d)', pair), 'Interpreter','latex','FontSize',14);
    legend('Location','northwest','Interpreter','latex');
    set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'XMinorGrid','on','YMinorGrid','on');
    set(gca, 'XScale', 'log');
    
    %% Right subplot: Outer Scaling
    subplot(1,2,2);
    hold on;
    for j = 1:2
        idx = 2*(pair-1) + j;
        plot(results(idx).y_over_delta, results(idx).U_over_Uinf, 'LineWidth', 2, ...
            'DisplayName', results(idx).baseString, 'Color', colors{j});
    end
    grid on;
    xlabel('$y/\delta$','Interpreter','latex','FontSize',12);
    ylabel('$U/U_{\infty}$','Interpreter','latex','FontSize',12);
    title(sprintf('Outer Scaling (Pair %d)', pair), 'Interpreter','latex','FontSize',14);
    legend('Location','best','Interpreter','latex');
    set(gca, 'FontSize', 12, 'LineWidth', 1.2);
    xlim([0,1])
    ylim([0,1])
    
    %% --- Add an annotation box with statistics for both cases ---
    % Prepare a stats string for the two cases in this pair.
% Construct the formatted string
statsStr = sprintf(['Case 1 (%s):\n' ...
    '$U_\\infty$ = %.2f m/s, $\\delta$ = %.0f mm\n' ...
    '$u_{\\tau}$ = %.2f m/s, $\\theta$ = %.4f m, $Re_{\\theta}$ = %.2e\n\n' ...
    'Case 2 (%s):\n' ...
    '$U_\\infty$ = %.2f m/s, $\\delta$ = %.0f mm\n' ...
    '$u_{\\tau}$ = %.2f m/s, $\\theta$ = %.4f m, $Re_{\\theta}$ = %.2e'], ...
    results(2*(pair-1)+1).baseString, results(2*(pair-1)+1).Uinf, results(2*(pair-1)+1).delta, ...
    results(2*(pair-1)+1).u_tau, results(2*(pair-1)+1).theta, results(2*(pair-1)+1).Re_theta, ...
    results(2*(pair-1)+2).baseString, results(2*(pair-1)+2).Uinf, results(2*(pair-1)+2).delta, ...
    results(2*(pair-1)+2).u_tau, results(2*(pair-1)+2).theta, results(2*(pair-1)+2).Re_theta);
    
    % Add an annotation textbox to the figure (position is in normalized figure units)
    annotation('textbox', [0.65, 0.15, 0.3, 0.3], 'String', statsStr, ...
        'Interpreter', 'latex', 'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 10);
end
