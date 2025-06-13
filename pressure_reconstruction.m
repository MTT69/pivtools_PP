function pressure_reconstruction(setup, CameraNo)
    % PRESSURE_RECONSTRUCTION - Reconstructs instantaneous pressure fields from PIV data using OMNIPOL
    % Inputs:
    %   setup - setup structure containing all parameters
    %   CameraNo - camera number
    
    % Water properties at room temperature (20°C)
    rho_water = 998.2;  % kg/m³
    nu_water = 1.004e-6; % m²/s (kinematic viscosity)
    
    if setup.pipeline.pressure_reconstruction
        for i = setup.instantaneous.runs
            disp(['Processing pressure reconstruction for run: ', num2str(i)]);
            
            % Define directories
            base_dir = setup.directory.base;
            data_dir = fullfile(base_dir, 'CalibratedPIV', num2str(setup.imProperties.imageCount), ...
                               ['Cam', num2str(CameraNo)], 'Instantaneous');
            stats_dir = fullfile(base_dir, 'Statistics', num2str(setup.imProperties.imageCount), ...
                                ['Cam' num2str(CameraNo)], 'Instantaneous', 'Calibrated');
            
            % Load mean statistics and coordinates
            window_size_str = [num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2))];
            MeanStats = load(fullfile(stats_dir, ['MeanStats' window_size_str '.mat']));
            Co_ords = load(fullfile(data_dir, 'Co_ords.mat'));
            
            % Get mean velocities and coordinates
            U_MEAN = MeanStats.mean_U;
            V_MEAN = MeanStats.mean_V;
            X_coords = Co_ords.Co_ords(i).x;
            
            % Calculate grid spacing (assuming uniform grid)
            dx = abs(X_coords(1,2) - X_coords(1,1));
            
            % Create output directory
            pressure_dir = fullfile(data_dir, 'Pressure');
            if ~exist(pressure_dir, 'dir')
                mkdir(pressure_dir);
            end
            
            % Process each image
            for imNo = 1:setup.imProperties.imageCount
                % Load current velocity data
                VelData = load(fullfile(data_dir, [num2str(sprintf(setup.instantaneous.nameConvention{1}, imNo))]));
                u = VelData.piv_result(i).ux;
                v = VelData.piv_result(i).uy;
                b_mask = VelData.piv_result(i).b_mask;
                
                % Create mask for OMNIPOL (ones with NaNs at masked locations)
                mask = ones(size(u));
                mask(b_mask) = NaN;
                mask(isnan(u) | isnan(v)) = NaN;
                
                % Apply NaN mask to velocities
                u(b_mask) = NaN;
                v(b_mask) = NaN;
                
                % Reconstruct pressure using OMNIPOL
                try
                    % Use 16 integration angles for better accuracy, with smoothing
                    P = OMNIPOL(u, v, U_MEAN, V_MEAN, dx, mask, nu_water, rho_water, 16, true, true);
                    
                    % Apply NaN mask to pressure field
                    P(b_mask) = NaN;
                    
                    % Save pressure field
                    pressure_result = struct();
                    pressure_result.P = P;
                    save_filename = fullfile(pressure_dir, [num2str(sprintf('%05d.mat', imNo))]);
                    save(save_filename, 'pressure_result');
                    
                    if mod(imNo, 100) == 0
                        disp(['Processed pressure for image ', num2str(imNo), ' of ', num2str(setup.imProperties.imageCount)]);
                    end
                    
                catch ME
                    disp(['Error processing image ', num2str(imNo), ': ', ME.message]);
                    continue;
                end
            end
            
            % Save reconstruction parameters
            pressure_params = struct();
            pressure_params.rho_water = rho_water;
            pressure_params.nu_water = nu_water;
            pressure_params.dx = dx;
            pressure_params.method = 'OMNIPOL';
            pressure_params.num_integration_angles = 16;
            pressure_params.smoothing_enabled = true;
            pressure_params.window_size = [setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)];
            pressure_params.total_images = setup.imProperties.imageCount;
            
            save(fullfile(pressure_dir, 'pressure_reconstruction_params.mat'), 'pressure_params');
            
            disp(['Pressure reconstruction completed for run ', num2str(i)]);
            disp(['Output saved to: ', pressure_dir]);
        end
    end
end
