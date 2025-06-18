function BL_analysis(setup, CameraNo, endpoint, x_loc, type, use_merged)
    if nargin < 5
        type = 'Instantaneous'; % Default for backward compatibility
    end
    if nargin < 6
        use_merged = false; % Default to traditional camera data
    end
    
    if setup.pipeline.batch_BL_variation 
        % Get data paths using centralized function
        paths = get_data_paths(setup, type, endpoint, CameraNo, use_merged);
        dataloc = paths.data_dir;
        statistics = paths.stats_dir;

        for i = setup.instantaneous.runs
            filename = fullfile(statistics, ['MeanStats' num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2)) '.mat']);
            meanData = load(filename);
            % Load coordinates
            Co_ords = load(fullfile(dataloc, 'Co_ords.mat'));
            x = Co_ords.Co_ords(i).x;
            y = Co_ords.Co_ords(i).y;
            
            % Initialize storage for batches
            num_batches = setup.imProperties.imageCount / setup.imProperties.caseImages;
            batch_U_profiles = cell(length(x_loc), num_batches);
            batch_u_prime_squared = cell(length(x_loc), num_batches);
            batch_v_prime_squared = cell(length(x_loc), num_batches);
            batch_u_prime_v_prime = cell(length(x_loc), num_batches);
            
            % Positive y indices
            positive_y_indices = y(:,1) > 0;
            y_positive = y(positive_y_indices, 1);
            
            % Loop over batches
            for b = 1:num_batches
                batch_U_sum = zeros(size(x));
                batch_u_prime_u_prime_sum = zeros(size(x));
                batch_v_prime_v_prime_sum = zeros(size(x));
                batch_u_prime_v_prime_sum = zeros(size(x));
                
                % Loop over images within the current batch
                for imNo = (b-1)*setup.imProperties.caseImages + 1 : b*setup.imProperties.caseImages
                    fileName = sprintf('%05d.mat', imNo);
                    VelData = load(fullfile(dataloc, fileName));
                    
                    % Accumulate velocity data
                    U = VelData.piv_result(i).ux;
                    V = VelData.piv_result(i).uy;
                    batch_U_sum = batch_U_sum + U;
                    
                    % Calculate fluctuations from overall mean
                    u_prime = U - meanData.mean_U;
                    v_prime = V - meanData.mean_V;
                    
                    % Accumulate turbulence quantities
                    batch_u_prime_u_prime_sum = batch_u_prime_u_prime_sum + u_prime.^2;
                    batch_v_prime_v_prime_sum = batch_v_prime_v_prime_sum + v_prime.^2;
                    batch_u_prime_v_prime_sum = batch_u_prime_v_prime_sum + u_prime.*v_prime;
                end
                
                % Compute batch averages
                batch_U_mean = batch_U_sum / setup.imProperties.caseImages;
                batch_u_prime_u_prime_mean = batch_u_prime_u_prime_sum / setup.imProperties.caseImages;
                batch_v_prime_v_prime_mean = batch_v_prime_v_prime_sum / setup.imProperties.caseImages;
                batch_u_prime_v_prime_mean = batch_u_prime_v_prime_sum / setup.imProperties.caseImages;
                
                % Calculate u_inf for this batch (99th percentile)
                batch_u_inf = prctile(batch_U_mean(positive_y_indices, :), 99, 'all');
                
                % Extract profiles at x_locations for this batch
                for x_idx = 1:length(x_loc)
                    [~, x_index] = min(abs(x(1,:) - x_loc(x_idx)));
                    
                    % Average 3 vectors either side (7 total)
                    col_range = max(1, x_index-3):min(size(x, 2), x_index+3);
                    
                    % Extract velocity profile
                    batch_U_profiles{x_idx, b} = mean(batch_U_mean(positive_y_indices, col_range), 2) / batch_u_inf;
                    
                    % Extract turbulence profiles
                    batch_u_prime_squared{x_idx, b} = mean(batch_u_prime_u_prime_mean(positive_y_indices, col_range), 2) / batch_u_inf^2;
                    batch_v_prime_squared{x_idx, b} = mean(batch_v_prime_v_prime_mean(positive_y_indices, col_range), 2) / batch_u_inf^2;
                    batch_u_prime_v_prime{x_idx, b} = mean(batch_u_prime_v_prime_mean(positive_y_indices, col_range), 2) / batch_u_inf^2;
                end
            end
            
            % Calculate overall average u_inf and delta
            overall_u_inf = prctile(meanData.mean_U(positive_y_indices, :), 99, 'all');
            [~, y_delta_indices] = min(abs(meanData.mean_U(positive_y_indices, :) - overall_u_inf), [], 1);
            y_delta_positions = y_positive(y_delta_indices);
            
            % Calculate momentum thickness for each batch and average
            batch_theta = zeros(length(x_loc), num_batches);
            avg_theta = zeros(length(x_loc), 1);
            
            for x_idx = 1:length(x_loc)
                [~, x_index] = min(abs(x(1,:) - x_loc(x_idx)));
                delta = y_delta_positions(x_index);
                
                % Calculate average profiles
                avg_U_profile = mean(cat(2, batch_U_profiles{x_idx, :}), 2);
                
                % Calculate momentum thickness for average profile
                [~, y_delta_idx] = min(abs(avg_U_profile - 1));
                y_restricted = flipud(y_positive(y_delta_idx:end));
                U_restricted = flipud(avg_U_profile(y_delta_idx:end));
                integrand = (U_restricted) .* (1 - U_restricted);
                avg_theta(x_idx) = trapz(y_restricted, integrand);
                
                % Calculate momentum thickness for each batch
                for b = 1:num_batches
                    U_batch = batch_U_profiles{x_idx, b};
                    [~, y_delta_idx] = min(abs(U_batch - 1));
                    if y_delta_idx < length(y_positive)
                        y_restricted = flipud(y_positive(y_delta_idx:end));
                        U_restricted = flipud(U_batch(y_delta_idx:end));
                        integrand = (U_restricted) .* (1 - U_restricted);
                        batch_theta(x_idx, b) = trapz(y_restricted, integrand);
                    end
                end
                
                % Create velocity profile plots for this x_location
                figure('Visible', 'off');
                set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
                hold on;
                
                colors = lines(num_batches + 1);
                
                % Plot individual batch profiles
                for b = 1:num_batches
                    plot(y_positive / delta, batch_U_profiles{x_idx, b}, ...
                        'Color', colors(b, :), 'LineWidth', 1, 'DisplayName', sprintf('Batch %d', b));
                end
                
                % Plot average profile
                plot(y_positive / delta, avg_U_profile, 'k-', 'LineWidth', 3, ...
                    'DisplayName', 'Average');
                
                xlabel('$y/\delta$', 'Interpreter', 'latex', 'FontSize', setup.figures.labelFontSize);
                ylabel('$U/U_\infty$', 'Interpreter', 'latex', 'FontSize', setup.figures.labelFontSize);
                title(sprintf('Velocity Profiles at x = %d mm', x_loc(x_idx)), ...
                      'Interpreter', 'latex', 'FontSize', setup.figures.titleFontSize);
                
                % Add annotations
                text(0.05, 0.95, sprintf('$U_\\infty = %.2f$ m/s', overall_u_inf), ...
                     'Units', 'normalized', 'Interpreter', 'latex', 'FontSize', setup.figures.legendFontSize);
                text(0.05, 0.90, sprintf('$\\delta = %.2f$ mm', delta), ...
                     'Units', 'normalized', 'Interpreter', 'latex', 'FontSize', setup.figures.legendFontSize);
                text(0.05, 0.85, sprintf('$\\theta_{avg} = %.2f$ mm', avg_theta(x_idx)), ...
                     'Units', 'normalized', 'Interpreter', 'latex', 'FontSize', setup.figures.legendFontSize);
                
                legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', setup.figures.legendFontSize);
                set(gca, 'FontSize', setup.figures.axisFontSize);
                grid on;
                xlim([0, 1]);
                ylim([0, 1]);
                hold off;
                
                % Save velocity profile plot
                vel_filename = fullfile(statistics, sprintf('BatchVelocityProfiles_x%d_run%d', x_loc(x_idx), i));
                saveas(gcf, [vel_filename, '.fig']);
                saveas(gcf, [vel_filename, '.png']);
                close(gcf);
                
                % Create Reynolds stress profile plots for this x_location
                figure('Visible', 'off');
                set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
                hold on;
                
                % Calculate average stress profiles
                avg_u_prime_squared = mean(cat(2, batch_u_prime_squared{x_idx, :}), 2);
                avg_v_prime_squared = mean(cat(2, batch_v_prime_squared{x_idx, :}), 2);
                avg_u_prime_v_prime = mean(cat(2, batch_u_prime_v_prime{x_idx, :}), 2);
                
                % Plot individual batch stress profiles
                for b = 1:num_batches
                    plot(y_positive / delta, batch_u_prime_squared{x_idx, b}, '--', ...
                        'Color', colors(b, :), 'LineWidth', 1, 'DisplayName', sprintf('$u''u''$ Batch %d', b));
                    plot(y_positive / delta, batch_v_prime_squared{x_idx, b}, '-.', ...
                        'Color', colors(b, :), 'LineWidth', 1, 'DisplayName', sprintf('$v''v''$ Batch %d', b));
                    plot(y_positive / delta, -batch_u_prime_v_prime{x_idx, b}, ':', ...
                        'Color', colors(b, :), 'LineWidth', 1, 'DisplayName', sprintf('$-u''v''$ Batch %d', b));
                end
                
                % Plot average stress profiles
                plot(y_positive / delta, avg_u_prime_squared, 'k--', 'LineWidth', 3, ...
                    'DisplayName', '$\overline{u''u''}$ Average');
                plot(y_positive / delta, avg_v_prime_squared, 'k-.', 'LineWidth', 3, ...
                    'DisplayName', '$\overline{v''v''}$ Average');
                plot(y_positive / delta, -avg_u_prime_v_prime, 'k:', 'LineWidth', 3, ...
                    'DisplayName', '$-\overline{u''v''}$ Average');
                
                xlabel('$y/\delta$', 'Interpreter', 'latex', 'FontSize', setup.figures.labelFontSize);
                ylabel('Reynolds Stress$/U_\infty^2$', 'Interpreter', 'latex', 'FontSize', setup.figures.labelFontSize);
                title(sprintf('Reynolds Stress Profiles at x = %d mm', x_loc(x_idx)), ...
                      'Interpreter', 'latex', 'FontSize', setup.figures.titleFontSize);
                
                % Add annotations
                text(0.05, 0.95, sprintf('$U_\\infty = %.2f$ m/s', overall_u_inf), ...
                     'Units', 'normalized', 'Interpreter', 'latex', 'FontSize', setup.figures.legendFontSize);
                text(0.05, 0.90, sprintf('$\\delta = %.2f$ mm', delta), ...
                     'Units', 'normalized', 'Interpreter', 'latex', 'FontSize', setup.figures.legendFontSize);
                text(0.05, 0.85, sprintf('$\\theta_{avg} = %.2f$ mm', avg_theta(x_idx)), ...
                     'Units', 'normalized', 'Interpreter', 'latex', 'FontSize', setup.figures.legendFontSize);
                
                legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', setup.figures.legendFontSize);
                set(gca, 'FontSize', setup.figures.axisFontSize);
                grid on;
                xlim([0, 1]);
                hold off;
                
                % Save Reynolds stress plot
                stress_filename = fullfile(statistics, sprintf('BatchReynoldsStress_x%d_run%d', x_loc(x_idx), i));
                saveas(gcf, [stress_filename, '.fig']);
                saveas(gcf, [stress_filename, '.png']);
                close(gcf);
            end
            
            % Save batch analysis data
            batch_data = struct();
            batch_data.x_locations = x_loc;
            batch_data.u_inf = overall_u_inf;
            batch_data.delta_positions = y_delta_positions;
            batch_data.avg_theta = avg_theta;
            batch_data.batch_theta = batch_theta;
            batch_data.batch_U_profiles = batch_U_profiles;
            batch_data.batch_stress_profiles.u_prime_squared = batch_u_prime_squared;
            batch_data.batch_stress_profiles.v_prime_squared = batch_v_prime_squared;
            batch_data.batch_stress_profiles.u_prime_v_prime = batch_u_prime_v_prime;
            
            batch_filename = fullfile(statistics, sprintf('BatchBLVariation_run%d.mat', i));
            save(batch_filename, '-struct', 'batch_data');
        end
    end
end