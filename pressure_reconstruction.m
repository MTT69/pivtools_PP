function pressure_reconstruction(setup, CameraNo)
    % PRESSURE_RECONSTRUCTION - Reconstructs instantaneous pressure fields from PIV data
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
        U_MEAN = MeanStats.U_mean;
        V_MEAN = MeanStats.V_mean;
        X = Co_ords.Co_ords(i).x;
        Y = Co_ords.Co_ords(i).y;
        
       
        % Time step
        dt = setup.imProperties.dt;
        
        % Create output directory
        pressure_dir = fullfile(data_dir, 'Pressure');
        if ~exist(pressure_dir, 'dir')
            mkdir(pressure_dir);
        end
        
        % Initialize storage for previous velocities (for time derivatives)
        U_prev = [];
        V_prev = [];
        U_curr = [];
        V_curr = [];
        
        % Process each image
        for imNo = 1:setup.imProperties.imageCount
            % Load current velocity data
            VelData = load(fullfile(data_dir, [num2str(sprintf(setup.instantaneous.nameConvention{1}, imNo))]));
            U_next = VelData.piv_result(i).ux;
            V_next = VelData.piv_result(i).uy;
            b_mask = VelData.piv_result(i).b_mask;
            
            % Apply NaN mask to velocities
            U_next(b_mask) = NaN;
            V_next(b_mask) = NaN;
            
            % Calculate time derivatives (need at least 2 time steps)
            if imNo == 1
                % First image: use forward difference for next calculation
                U_prev = U_next;
                V_prev = V_next;
                continue;
            elseif imNo == 2
                % Second image: use forward difference
                U_curr = U_prev;
                V_curr = V_prev;
                dUdt = (U_next - U_curr) / dt;
                dVdt = (V_next - V_curr) / dt;
            else
                % Use central difference for interior points
                dUdt = (U_next - U_prev) / (2 * dt);
                dVdt = (V_next - V_prev) / (2 * dt);
            end
            
            % Apply NaN mask to time derivatives
            dUdt(b_mask) = NaN;
            dVdt(b_mask) = NaN;
            
            % Reconstruct pressure using EU_FDM
            try
                P = EU_FDM(X, Y, U_curr, V_curr, U_MEAN, V_MEAN, dUdt, dVdt, rho_water, nu_water);
                
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
            
            % Update previous velocities for next iteration
            if imNo >= 2
                U_prev = U_curr;
                V_prev = V_curr;
                U_curr = U_next;
                V_curr = V_next;
            end
        end
        
        % Save reconstruction parameters
        pressure_params = struct();
        pressure_params.rho_water = rho_water;
        pressure_params.nu_water = nu_water;
        pressure_params.dt = dt;
        pressure_params.method = 'EU_FDM';
        pressure_params.window_size = [setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)];
        pressure_params.total_images = setup.imProperties.imageCount;
        
        save(fullfile(pressure_dir, 'pressure_reconstruction_params.mat'), 'pressure_params');
        
        disp(['Pressure reconstruction completed for run ', num2str(i)]);
        disp(['Output saved to: ', pressure_dir]);
    end
    end
end
