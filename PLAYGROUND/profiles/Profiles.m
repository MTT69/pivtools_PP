close all
clear variables
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

scalefactor = 9.53;
N=18000;
profiles_all=struct();

for base_idx = 1:length(base_dir)
    Base = base_dir{base_idx};
    % Extract the base string after the last backslash
    [~, baseString, ext] = fileparts(Base);

    % Extract dt from the base directory name
    dt_str = regexp(baseString, '\d+dt', 'match');
    dt = str2double(dt_str{1}(1:end-2)) * 10^(-6);

    % Print the file name and dt to the terminal
    fprintf('Processing file: %s with dt: %.6f\n', [baseString], dt);
    Camera = 1;
    window = [16 16];
    i = 5;
    lower_limit = -0.25;
    upper_limit = 0.25;

    calibrated = fullfile(Base, 'CalibratedPIV', num2str(N), ['Cam' num2str(Camera)], 'Instantaneous');
    Co_ords = load(fullfile(calibrated, "Co_ords.mat"));
    VelData = load(fullfile(calibrated, "00001.mat"));
    statistics = fullfile(Base, 'Statistics', num2str(N), ['Cam' num2str(Camera)], 'Instantaneous', 'Calibrated');
    filename = fullfile(statistics, ['MeanStats' num2str(window(1)) 'x' num2str(window(2)) '.mat']);
    meanData = load(filename);
    ycorners = [Co_ords.Co_ords(i).y(1,1), Co_ords.Co_ords(i).y(end, end)];
    xcorners = [Co_ords.Co_ords(i).x(1,1), Co_ords.Co_ords(i).x(end, end)];
    b_mask = VelData.piv_result(i).b_mask;

    x_locations = [-100, -5, 50];

    % Extract the column vector of velocities at specified x locations and filter positive y values

    velocity_profiles = cell(length(x_locations), 1);
    positive_y_indices = Co_ords.Co_ords(i).y(:,1) > 0;
    percentile_speed = prctile(meanData.mean_U(positive_y_indices, :), 99, 'all');
    disp(["U-inf", percentile_speed])
    x_index= cell(1,3);
    for idx = 1:length(x_locations)
        x_loc = x_locations(idx);
        [~, x_index{idx}] = min(abs(Co_ords.Co_ords(i).x(1,:) - x_loc));
        
        % Define the range of columns to average
        col_range = max(1, x_index{idx}-2):min(size(Co_ords.Co_ords(i).x, 2), x_index{idx}+2);
        
        % Extract and normalize velocity data at the current x location range
        velocity_profiles{idx} = mean(meanData.mean_U(positive_y_indices, col_range), 2) / percentile_speed;
    end

    % Calculate the 99th percentile speed for all mean U values at locations where y > 0

    % Find the y position at which the percentile speed is achieved
    [~, y_percentile_index] = min(abs(meanData.mean_U(positive_y_indices, :) - percentile_speed), [], 1);
    y_percentile_position = Co_ords.Co_ords(i).y(positive_y_indices, 1);
    y_percentile_position = y_percentile_position(y_percentile_index);

    fprintf('delta: %.4f mm\n', y_percentile_position(x_index{2}));
    fprintf('U_inf: %.4f m/s\n', percentile_speed);

    % Calculate the momentum thickness (theta) for each velocity profile
    momentum_thickness = zeros(length(x_locations), 1);
    delta_star = zeros(length(x_locations), 1);

    for idx = 1:length(x_locations)
        U = velocity_profiles{idx}; % Velocity profile at the current x location
        y = Co_ords.Co_ords(i).y(positive_y_indices, 1); % Corresponding y coordinates
        U_inf = 1; % Freestream velocity (99th percentile speed)
        
        % Find the index of the closest velocity to U_inf
        [~, y_delta_index] = min(abs(U - U_inf));
        y_delta = y(y_delta_index); % y-coordinate corresponding to U_inf

        % Find the index of the y-coordinate closest to zero
        [~, y_zero_index] = min(abs(y));
        
        % Restrict y and U to the range [y_zero, y_delta]
        y_restricted = flipud(y(y_delta_index: y_zero_index));
        U_restricted = flipud(U(y_delta_index: y_zero_index));
        
        % Calculate the integrand for momentum thickness
        integrand = (U_restricted ./ U_inf) .* (1 - (U_restricted ./ U_inf));
        delta_integrand = (1 - (U_restricted ./ U_inf));
        
        % Perform numerical integration using the trapezoidal rule
        theta = trapz(y_restricted, integrand);
        delta_star(idx) = trapz(y_restricted,delta_integrand);
        momentum_thickness(idx) = theta;
    end
     


    % Example of printing out the momentum thickness and displacement thickness
    fprintf('theta = @x = %d mm -->  %.4f mm \n', x_locations(2), momentum_thickness(2));
    fprintf('delta* = @x = %d mm -->  %.4f mm \n', x_locations(2), delta_star(2));
    % Define the range for fitting
    y_min = 0.05 * y_percentile_position(x_index{2});
    y_max = 0.15 * y_percentile_position(x_index{2});

    % Filter the data within the specified range
    fit_indices = (y >= y_min) & (y <= y_max);
    y_fit = y(fit_indices)*10^-3; %convert to m
    U_fit = meanData.mean_U(fit_indices,x_index{2});

    % Define the von Kármán constant and empirical constant
    kappa = 0.39;
    B = 4.3;
    nu = 1*10^-6;

    u_tau_initial_guess = 0.05* percentile_speed;

    % Define the function to minimize
    log_law = @(u_tau) ((u_tau / kappa) * log(y_fit * u_tau / nu) + B*u_tau - U_fit);

    % Use least squares solver to find u_tau with Levenberg-Marquardt algorithm, lower tolerance, and increased iterations
    lsqoptions = optimoptions('lsqnonlin', 'Display', 'off', 'Algorithm', 'levenberg-marquardt', 'TolFun', 1e-10, 'TolX', 1e-10, 'MaxIter', 1000, 'MaxFunEvals', 2000);
    [u_tau, resnorm, residual, exitflag, output] = lsqnonlin(log_law, u_tau_initial_guess, [], [], lsqoptions);

    fprintf('Estimated u_tau: %.4f\n', u_tau);
    fprintf('Residual norm: %.4e\n', resnorm);

    % Calculate the log law profile using the estimated u_tau
    y_log = linspace(min(y_fit), max(y_fit), 100);
    U_log = (u_tau / kappa) * log(y_log * u_tau / nu) + B * u_tau;

    

    % Calculate and print the momentum thickness Reynolds number (Re_theta)
    nu = 1e-6; % Kinematic viscosity (example value, adjust as needed)
    Re_theta = percentile_speed * momentum_thickness(2)*10^-3 / nu;

    fprintf('Re_theta @ x = %d --> %.4e\n', x_locations(2), Re_theta);

    % Mark the x = 0 position momentum thickness on the graph
    figure(1);
    set(gcf, 'Visible', 'off');
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    hold on;
    colors = lines(length(x_locations)); % Generate distinct colors for plots
    momentum_thickness_line = []; % Placeholder for the momentum thickness line handle

    for idx = 1:length(x_locations)
        % Plot velocity profiles with markers
        plot(Co_ords.Co_ords(i).y(positive_y_indices, 1)/y_percentile_position(x_index{2}), velocity_profiles{idx}, ...
            'LineWidth', 2, 'Color', colors(idx, :), 'Marker', 'o', ...
            'DisplayName', ['x = ' num2str(x_locations(idx))]);
    end

    % Plot the estimated log law profile
    plot(y_log*10^3/y_percentile_position(x_index{2}) , U_log/percentile_speed , 'k--', 'LineWidth', 2, 'DisplayName', 'Log Law Fit');
    ylim([0,1])
    xlim([0,1])
    hold off;
    xlabel('$Y/\delta$ ', 'Interpreter', 'latex');
    ylabel('Velocity ($U/U_\infty$)', 'Interpreter', 'latex');
%     set(gca, 'XScale', 'log', 'TickLabelInterpreter', 'latex');

    % Add legend with the correct entries for velocity profiles and momentum thickness
    legend('show', 'Location', 'Best', 'Interpreter', 'latex'); 
    title('Velocity Profiles at Specified X Locations', 'Interpreter', 'latex');
    set(gca, 'FontSize', 16);
    grid on;

    y_plus = y*(10^-3) * u_tau / nu;
    u_plus = velocity_profiles{2}*percentile_speed / u_tau;

    % Create a new figure for the u+ vs y+ plot
    figure(6);
    set(gcf, 'Visible', 'off');
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    hold on;

    % Plot u+ against y+
    plot(y_plus, u_plus, 'bo-', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', '$u^+$');

    % % Plot the log law profile
    % y_plus_log = linspace(min(y_plus), max(y_plus), 100);
    % u_plus_log = (1 / kappa) * log(y_plus_log) + B;
    % plot(y_plus_log, u_plus_log, 'k--', 'LineWidth', 2, 'DisplayName', 'Log Law Fit');

    hold off;
    xlabel('$y^+$', 'Interpreter', 'latex');
    ylabel('$u^+$', 'Interpreter', 'latex');
    title('$u^+$ vs $y^+$', 'Interpreter', 'latex');
    legend('show', 'Location', 'Best', 'Interpreter', 'latex');
    set(gca, 'XScale', 'log', 'FontSize', 16, 'TickLabelInterpreter', 'latex');
    grid on;
    ylim([0,25])

    u_prime_squared_profiles = cell(length(x_locations), 1);
    v_prime_squared_profiles = cell(length(x_locations), 1);
    u_prime_v_prime_profiles = cell(length(x_locations), 1);
    positive_y_indices = Co_ords.Co_ords(i).y(:,1) > 0;

    for idx = 1:length(x_locations)
        x_loc = x_locations(idx);
        [~, x_loc] = min(abs(Co_ords.Co_ords(i).x(1,:) - x_loc));
        
        % Define the range of columns to average
        col_range = max(1, x_loc-2):min(size(Co_ords.Co_ords(i).x, 2), x_loc+2);
        
        % Extract and normalize turbulence data at the current x location range
        u_prime_squared_profiles{idx} = mean(meanData.U_prime_Uprime_mean(positive_y_indices, col_range), 2) / percentile_speed^2;
        v_prime_squared_profiles{idx} = mean(meanData.V_prime_Vprime_mean(positive_y_indices, col_range), 2) / percentile_speed^2;
        u_prime_v_prime_profiles{idx} = mean(meanData.U_prime_Vprime_mean(positive_y_indices, col_range), 2) / percentile_speed^2;
    end

    y_positive = Co_ords.Co_ords(i).y(positive_y_indices, 1); % y values where y > 0

    % Plot all turbulence profiles on the same graph
    figure;
    set(gcf, 'Visible', 'off');
    hold on;
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

    % Define line styles and markers for different stresses
    line_styles = {'-', '--', '-.'}; % Solid for u'^2, dashed for v'^2, dash-dot for u'v'
    markers = {'o', 's', 'd'}; % Circle for u'^2, square for v'^2, diamond for u'v'

    % Plot each stress profile for each x location
    for idx = 1:length(x_locations)
        % Plot u'^2 with a solid line
        plot(y_positive, u_prime_squared_profiles{idx}, ...
            'LineStyle', line_styles{1}, 'Color', colors(idx, :), 'LineWidth', 1, ...
            'Marker', markers{1}, 'DisplayName', sprintf('$u''^2$, x = %d', x_locations(idx)));
        
        % Plot v'^2 with a dashed line
        plot(y_positive, v_prime_squared_profiles{idx}, ...
            'LineStyle', line_styles{2}, 'Color', colors(idx, :), 'LineWidth', 1, ...
            'Marker', markers{2}, 'DisplayName', sprintf('$v''^2$, x = %d', x_locations(idx)));
        
        % Plot u'v' with a dash-dot line
        plot(y_positive,  -10*u_prime_v_prime_profiles{idx}, ...
            'LineStyle', line_styles{3}, 'Color', colors(idx, :), 'LineWidth', 1, ...
            'Marker', markers{3}, 'DisplayName', sprintf('$-10u''v''$, x = %d', x_locations(idx)));
    end

    % Customize the graph
    xlabel('Y (mm)', 'Interpreter', 'latex');
    ylabel('Stress / $U_\infty^2$', 'Interpreter', 'latex');
    title('Normalised Turbulence Stress Profiles', 'Interpreter', 'latex');
    legend('show', 'Location', 'Best', 'Interpreter', 'latex'); % Display a legend for all lines
    grid on;
    set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
    hold off;
    ylim([-2*10^-3, 1.4*10^-2])

    figure(4);
    set(gcf, 'Visible', 'off');
    mean_U = meanData.mean_U;
    mean_V = meanData.mean_V;
    normalised_vorticity = meanData.mean_Vorticity*momentum_thickness(2)/percentile_speed;
    center_col = floor(size(mean_U, 2) / 2);
    half_width = floor(size(mean_U, 1) / 2);
    U_cropped = mean_U(:, center_col - half_width : center_col + half_width);
    V_cropped = mean_V(:, center_col - half_width : center_col + half_width);
    x_co_ords_cropped = Co_ords.Co_ords(i).x(:, center_col - half_width : center_col + half_width);
    y_co_ords_cropped = Co_ords.Co_ords(i).y(:, center_col - half_width : center_col + half_width);
    xcorners_cropped = [x_co_ords_cropped(1,1), x_co_ords_cropped(1,end)];
    ycorners_cropped = [y_co_ords_cropped(1,1), y_co_ords_cropped(end, end)];
    b_mask_cropped = b_mask(:, center_col - half_width : center_col + half_width);
    vorticity_cropped = normalised_vorticity(:, center_col - half_width : center_col + half_width);
    upsampled_mean_u = imresize(U_cropped, 5, "bicubic");
    upsampled_mean_v = imresize(V_cropped, 5, "bicubic");
    upsampled_b_mask = imresize(b_mask_cropped, 5, "nearest");
    upsampled_vorticity = imresize(vorticity_cropped, 5, "bicubic");

    upsample_cat_v = cat(3, -upsampled_mean_v, upsampled_mean_u);
    upsample_cat_v = perform_vf_normalization(upsample_cat_v);

    options.spot_size = 2;
    options.flow_correction = 1;

    lic = perform_lic(upsample_cat_v, 12, options);
    lic(upsampled_b_mask==1) =0;
    upsampled_vorticity(upsampled_b_mask==1)=0;

    ax1 = axes;
    A = imagesc(xcorners_cropped, ycorners_cropped,lic);
    colormap(ax1, gray);
    daspect(ax1, [1 1 1]);
    alphaData_lic = lic ~= 0; % True for non-zero values, false for zeros
    set(A, 'AlphaData', alphaData_lic);
    set(gca, 'YDir', 'normal', 'FontSize', 16, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex'); % Thicker axes, LaTeX labels            colormap(ax2, [0 0 0]);
    hold on;
    ax2 = axes;
    B = imagesc(xcorners_cropped, ycorners_cropped,upsampled_vorticity,[lower_limit, upper_limit]);
    colormap(ax2, redbluezero(lower_limit, upper_limit));
    daspect(ax2, [1 1 1]);
    colorbar;
    % add a colorbar description string of normalised vorticity (omega*theta/U_inf)
    c = colorbar;
    c.Label.String = 'Normalised Vorticity ($\omega\theta/U_\infty$)';
    %make the c.label have a latex interpreter
    c.Label.Interpreter = 'latex';
    alphaData = ~isnan(upsampled_vorticity);
    set(B, 'AlphaData', alphaData * 0.5); % Use 0.5 transparency for non-zero values
    set(gca, 'YDir', 'normal', 'FontSize', 16, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex'); % Thicker axes, LaTeX labels            colormap(ax2, [0 0 0]);

    % % Align axes
    linkaxes([ax1, ax2]);
    % 
    % % Hide the second set of axes
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    ax2.HitTest = 'off'; % Ignore mouse clicks on ax2
    ax2.PickableParts = 'none'; % Prevent ax2 from capturing interactions

    ax2.Position = ax1.Position;
    daspect([1 1 1]);
    % 
    % % Remove extra space between plots


    ax1.XLim = [-25, 45];
    ax2.XLim = [-25, 45];
    ax1.YLim = [-50, 10];
    ax2.YLim = [-50, 10];
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);


    % Calculate the scaled mean_U

    scaled_mean_U = meanData.mean_U * (scalefactor/10.^-3*dt);

    % Create a new figure for the scaled mean_U and signal-to-noise ratio
    figure(5);
    set(gcf, 'Visible', 'off');
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

    % Plot the scaled mean_U using the redbluezero colormap
    subplot(1, 2, 1);
    imagesc(xcorners, ycorners, abs(scaled_mean_U), [0,5]);
    colormap(parula);
    colorbar;
    %add a description to colorbar of scaled mean U in pixels
    c = colorbar;
    c.Label.String = 'abs Mean U (pixels)';
    c.Label.Interpreter = 'latex';
    title('Scaled Mean U', 'Interpreter', 'latex');
    xlabel('X (mm)', 'Interpreter', 'latex');
    ylabel('Y (mm)', 'Interpreter', 'latex');
    set(gca, 'YDir', 'normal', 'FontSize', 16, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex'); % Thicker axes, LaTeX labels  
    ax1 = gca;
    ax1.XLim = [-25, 45];
    ax1.YLim = [-50, 10];

    % Calculate the signal-to-noise ratio
    signal_to_noise = abs(scaled_mean_U) / 0.1;

    % Plot the signal-to-noise ratio
    subplot(1, 2, 2);
    imagesc(xcorners, ycorners, signal_to_noise, [0,15]);
    colormap("parula");
    colorbar;
    %add a description to colorbar of signal to noise ratio
    c = colorbar;
    c.Label.String = 'Signal-to-Noise Ratio';
    title('Signal-to-Noise Ratio', 'Interpreter', 'latex');
    xlabel('X (mm)', 'Interpreter', 'latex');
    ylabel('Y (mm)', 'Interpreter', 'latex');
    set(gca, 'YDir', 'normal', 'FontSize', 16, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex'); % Thicker axes, LaTeX labels  
    ax2 = gca;
    ax2.XLim = [-25, 45];
    ax2.YLim = [-50, 10];
    daspect([1 1 1])

    % 
    % 
    % % Save the velocity profiles figure
    figure(1)
    set(gcf, 'Visible', 'off');
    velocity_profiles_filename = fullfile(statistics, sprintf('VelocityProfiles_%dx%d', window(1), window(2)));
    saveas(gcf, join([velocity_profiles_filename, '.fig'], ''));
    saveas(gcf, join([velocity_profiles_filename, '.jpeg'], ''));
    % 
    % % Save the turbulence stress profiles figure
    figure(2); % Assuming the turbulence stress profiles figure is the second figure
    set(gcf, 'Visible', 'off');
    turbulence_profiles_filename = fullfile(statistics, sprintf('TurbulenceProfiles_%dx%d', window(1), window(2)));
    saveas(gcf, join([turbulence_profiles_filename, '.fig'], ''));
    saveas(gcf, join([turbulence_profiles_filename, '.jpeg'], ''));
    % 

    %save figure 4 with the name vorticity and the windowsizes as above
    figure(4)
    set(gcf, 'Visible', 'off');
    vorticity_filename = fullfile(statistics, sprintf('Vorticity_%dx%d', window(1), window(2)));
    saveas(gcf, join([vorticity_filename, '.fig'], ''));
    saveas(gcf, join([vorticity_filename, '.jpeg'], ''));

    %save figure 5 called SNR and windowsizes as above
    figure(5)
    set(gcf, 'Visible', 'off');
    SNR_filename = fullfile(statistics, sprintf('SNR_%dx%d', window(1), window(2)));
    saveas(gcf, join([SNR_filename, '.fig'], ''));
    saveas(gcf, join([SNR_filename, '.jpeg'], ''));

    %save figure 6 called utau_profiles and windowsizes as above
    figure(6)
    set(gcf, 'Visible', 'off');
    utau_profiles_filename = fullfile(statistics, sprintf('utau_profiles_%dx%d', window(1), window(2)));
    saveas(gcf, join([utau_profiles_filename, '.fig'], ''));
    saveas(gcf, join([utau_profiles_filename, '.jpeg'], ''));


    profile_data = struct();
    profile_data.u_inf = percentile_speed;
    profile_data.delta = y_percentile_position(x_index{2});
    profile_data.delta_star = delta_star(2);
    profile_data.u_tau = u_tau;
    profile_data.theta = momentum_thickness(2);

    profile_filename = fullfile(statistics, sprintf('profile_%dx%d.mat', window(1), window(2)));
    save(profile_filename, '-struct', 'profile_data');
        profile_data = struct();

    profiles_all(base_idx).u_inf = percentile_speed;
    profiles_all(base_idx).delta = y_percentile_position(x_index{2});
    profiles_all(base_idx).delta_star = delta_star(2);
    profiles_all(base_idx).u_tau = u_tau;
    profiles_all(base_idx).theta = momentum_thickness(2);

    close all
end
% Extract indices for 100 Hz and 250 Hz cases
%% Extract Shortened Legends
short_legends = cell(size(base_dir));
for i = 1:length(base_dir)
    tokens = regexp(base_dir{i}, '(\d+degree)(_?)(\d+light)(_?)(\d+hz)', 'tokens');
    if ~isempty(tokens)
        short_legends{i} = sprintf('%s', tokens{1}{1});
    else
        short_legends{i} = sprintf('Case %d', i); % Fallback name
    end
end

%% Separate Data by Frequency
delta_star_100Hz = [];
theta_100Hz = [];
legends_100Hz = {};

delta_star_250Hz = [];
theta_250Hz = [];
legends_250Hz = {};

for i = 1:length(base_dir)
    if contains(base_dir{i}, '100hz')
        delta_star_100Hz = [delta_star_100Hz, profiles_all(i).delta_star];
        theta_100Hz = [theta_100Hz, profiles_all(i).theta];
        legends_100Hz{end+1} = short_legends{i};
    elseif contains(base_dir{i}, '250hz')
        delta_star_250Hz = [delta_star_250Hz, profiles_all(i).delta_star];
        theta_250Hz = [theta_250Hz, profiles_all(i).theta];
        legends_250Hz{end+1} = short_legends{i};
    end
end

%% Plot Bar Chart for 100Hz Cases
figure;
bar([delta_star_100Hz; theta_100Hz]', 'grouped');
grid on;
set(gca, 'FontSize', 14, 'LineWidth', 1.5);
ylabel('mm', 'FontSize', 16);
xticklabels(legends_100Hz);
xtickangle(45);
legend({'\delta*', '\theta'}, 'FontSize', 14, 'Location', 'Best');
title('Boundary Layer Thickness (100Hz Cases)', 'FontSize', 16);

%% Plot Bar Chart for 250Hz Cases
figure;
bar([delta_star_250Hz; theta_250Hz]', 'grouped');
grid on;
set(gca, 'FontSize', 14, 'LineWidth', 1.5);
ylabel('mm', 'FontSize', 16);
xticklabels(legends_250Hz);
xtickangle(45);
legend({'\delta*', '\theta'}, 'FontSize', 14, 'Location', 'Best');
title('Boundary Layer Thickness (250Hz Cases)', 'FontSize', 16);

