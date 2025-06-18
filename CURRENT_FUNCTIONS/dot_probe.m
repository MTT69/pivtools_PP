function dot_probe(setup, CameraNo, endpoint, probe_x, probe_y, type, use_merged)
    % DOT_PROBE - Spatial averaging of velocity components over 3x3 pixel probe
    % Inputs:
    %   setup - setup structure
    %   CameraNo - camera number
    %   endpoint - endpoint string
    %   probe_x - x-coordinate of probe center in physical units
    %   probe_y - y-coordinate of probe center in physical units
    %   type - 'Instantaneous' or 'Ensemble'
    %   use_merged - logical, true to use merged data, false for traditional camera data
    
    if nargin < 6
        type = 'Instantaneous'; % Default for backward compatibility
    end
    if nargin < 7
        use_merged = false; % Default to traditional camera data
    end
    
    if setup.pipeline.dot_probe
        % Get data paths using centralized function
        paths = get_data_paths(setup, type, endpoint, CameraNo, use_merged);
        dataloc = paths.data_dir;
        
        fprintf('Performing dot probe analysis at physical coordinates (%.2f, %.2f) at %s\n', probe_x, probe_y, datetime('now'));
        
        for i = setup.instantaneous.runs
            VelData = load(fullfile(dataloc, sprintf(setup.instantaneous.nameConvention{1}, 1)));
            Co_ords = load(fullfile(dataloc, 'Co_ords.mat'));
            
            h = size(Co_ords.Co_ords(i).y, 1);
            w = size(Co_ords.Co_ords(i).x, 2);
            
            % Convert physical coordinates to pixel indices
            x_coords = Co_ords.Co_ords(i).x;
            y_coords = Co_ords.Co_ords(i).y;
            
            % Find closest pixel to physical coordinates
            [~, probe_x_idx] = min(abs(x_coords(1,:) - probe_x));
            [~, probe_y_idx] = min(abs(y_coords(:,1) - probe_y));
            
            fprintf('Physical coordinates (%.2f, %.2f) correspond to pixel indices (%d, %d)\n', ...
                    probe_x, probe_y, probe_x_idx, probe_y_idx);
            
            % Create 3x3 probe mask (center pixel + 8 surrounding pixels)
            probe_mask = false(h, w);
            for dy = -1:1
                for dx = -1:1
                    y_idx = probe_y_idx + dy;
                    x_idx = probe_x_idx + dx;
                    if y_idx >= 1 && y_idx <= h && x_idx >= 1 && x_idx <= w
                        probe_mask(y_idx, x_idx) = true;
                    end
                end
            end
            
            % Initialize time history storage
            num_batches = setup.imProperties.imageCount / setup.imProperties.caseImages;
            u_time_history = cell(num_batches, 1);
            v_time_history = cell(num_batches, 1);
            
            for batch = 1:(setup.imProperties.imageCount/setup.imProperties.caseImages)
                fprintf('Processing batch %d of %d for pass %d at %s\n', batch, num_batches, i, datetime('now'));
                
                % Preallocate for current batch
                u_batch = zeros(setup.imProperties.caseImages, 1);
                v_batch = zeros(setup.imProperties.caseImages, 1);
                
                % Loop over images in batch
                parfor idx = 1:setup.imProperties.caseImages
                    ImNo = (batch-1)*setup.imProperties.caseImages + idx;
                    
                    % Load velocity data
                    VelData = load(fullfile(dataloc, sprintf(setup.instantaneous.nameConvention{1}, ImNo)));
                    u = VelData.piv_result(i).ux;
                    v = VelData.piv_result(i).uy;
                    
                    % Spatially average over 3x3 probe region
                    u_probe_values = u(probe_mask);
                    v_probe_values = v(probe_mask);
                    
                    % Remove NaN values before averaging
                    u_probe_values = u_probe_values(~isnan(u_probe_values));
                    v_probe_values = v_probe_values(~isnan(v_probe_values));
                    
                    u_batch(idx) = mean(u_probe_values);
                    v_batch(idx) = mean(v_probe_values);
                end
                
                % Store batch time history
                u_time_history{batch} = u_batch;
                v_time_history{batch} = v_batch;
            end
            
            % Save results
            if use_merged
                output_dir = fullfile(setup.directory.base, 'Statistics', 'Merged', type, endpoint, 'DotProbe');
            else
                output_dir = fullfile(setup.directory.base, 'Statistics', num2str(setup.imProperties.imageCount), ...
                                     ['Cam' num2str(CameraNo)], type, endpoint, 'DotProbe');
            end
            if ~exist(output_dir, 'dir')
                mkdir(output_dir);
            end
            
            filename = sprintf('DotProbe_x%.2f_y%.2f_pass%d.mat', probe_x, probe_y, i);
            save(fullfile(output_dir, filename), 'u_time_history', 'v_time_history', ...
                 'probe_x', 'probe_y', 'probe_x_idx', 'probe_y_idx', 'num_batches', '-v7.3');
            
            % Create vertical subplots
            figure('Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
            
            % U component subplot
            subplot(2, 1, 1);
            hold on;
            colors = lines(num_batches);
            for batch = 1:num_batches
                time_vector = (1:length(u_time_history{batch})) + (batch-1)*setup.imProperties.caseImages;
                plot(time_vector, u_time_history{batch}, 'Color', colors(batch,:), 'LineWidth', 1.5);
            end
            ylabel('U Velocity', 'FontSize', setup.figures.labelFontSize);
            title(sprintf('U Component Time History - Probe at (%.2f, %.2f)', probe_x, probe_y), ...
                  'FontSize', setup.figures.titleFontSize);
            grid on; box on;
            set(gca, 'FontSize', setup.figures.axisFontSize);
            legend(arrayfun(@(x) sprintf('Batch %d', x), 1:num_batches, 'UniformOutput', false), ...
                   'Location', 'best');
            
            % V component subplot
            subplot(2, 1, 2);
            hold on;
            for batch = 1:num_batches
                time_vector = (1:length(v_time_history{batch})) + (batch-1)*setup.imProperties.caseImages;
                plot(time_vector, v_time_history{batch}, 'Color', colors(batch,:), 'LineWidth', 1.5);
            end
            xlabel('Image Number', 'FontSize', setup.figures.labelFontSize);
            ylabel('V Velocity', 'FontSize', setup.figures.labelFontSize);
            title(sprintf('V Component Time History - Probe at (%.2f, %.2f)', probe_x, probe_y), ...
                  'FontSize', setup.figures.titleFontSize);
            grid on; box on;
            set(gca, 'FontSize', setup.figures.axisFontSize);
            legend(arrayfun(@(x) sprintf('Batch %d', x), 1:num_batches, 'UniformOutput', false), ...
                   'Location', 'best');
            
            % Save figure
            fig_name = sprintf('DotProbe_TimeHistory_x%.2f_y%.2f_pass%d', probe_x, probe_y, i);
            saveas(gcf, fullfile(output_dir, [fig_name, '.jpg']));
            saveas(gcf, fullfile(output_dir, [fig_name, '.epsc']));
            saveas(gcf, fullfile(output_dir, [fig_name, '.fig']));
            close(gcf);
            
            fprintf('Dot probe analysis completed for pass %d at %s\n', i, datetime('now'));
        end
    end
end
