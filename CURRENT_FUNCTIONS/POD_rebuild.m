function POD_rebuild(setup, CameraNo, domain, type, use_merged)
    % POD_rebuild - Reconstructs flow field using POD modes based on noise analysis
    % Inputs:
    %   setup - setup structure containing all parameters
    %   CameraNo - camera number
    %   domain - 'Above', 'Below', or 'Full'
    %   type - 'Instantaneous', 'Calibrated', or 'Merged'
    %   use_merged - logical, true to use merged data, false for traditional camera data
    
    if nargin < 3
        domain = 'Below'; % Default to cavity-only reconstruction
    end
    if nargin < 4
        type = 'Instantaneous'; % Default for backward compatibility
    end
    if nargin < 5
        use_merged = false; % Default to traditional camera data
    end
    
    if setup.pipeline.POD_rebuild
       
    for i = setup.instantaneous.runs
        disp(['Processing run: ', num2str(i)]);
        
        % Get data paths using centralized function
        paths = get_data_paths(setup, type, '', CameraNo, use_merged);
        stats_dir = paths.stats_dir;
        data_dir = paths.data_dir;
        
        % Load mean statistics
        window_size_str = [num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2))];
        MeanStats = load(fullfile(stats_dir, 'meanStats', ['MeanStats' window_size_str '.mat']));
        
        % Load coordinate data
        Co_ords = load(fullfile(data_dir, 'Co_ords.mat'));
        VelData = load(fullfile(data_dir, [num2str(sprintf(setup.instantaneous.nameConvention{1}, 1))]));
        b_mask = VelData.piv_result(i).b_mask;
        
        % Get coordinates
        y = Co_ords.Co_ords(i).y;
        
        % --- Domain logic copied from POST_POD_multi ---
        mask = b_mask;
        row_average = mean(mask, 2);
        index = find(row_average > 0.8, 1);

        if strcmp(domain, 'Above')
            mask(index:end, :) = 1;
        elseif strcmp(domain, 'Below')
            mask(1:index-1, :) = 1;
        elseif strcmp(domain, 'Full')
            mask = mask;
        end
        % ------------------------------------------------

        % Find cavity region (y <= 0, not in mask)
        cavity_region = (y <= 0) & ~mask;
        
        % Calculate uncertainty in cavity region
        u_prime_squared_calibrated = MeanStats.U_prime_Uprime_mean;
        v_prime_squared_calibrated = MeanStats.V_prime_Vprime_mean;
        TKE_calibrated = 0.5 * (u_prime_squared_calibrated + v_prime_squared_calibrated);
        
        TKE_physical = TKE_calibrated;
        
        % Calculate total TKE in cavity region
        TKE_cavity = TKE_physical(cavity_region);
        Total_TKE = sum(TKE_cavity(:));
        
        % Calculate noise contribution
        uncertainty = 0.1 / (setup.imProperties.scaleFactor / 1e-3 * setup.imProperties.dt);
        cavity_area = sum(cavity_region(:));
        Noise_TKE = cavity_area * uncertainty^2;
        
        % Calculate SNR and signal percentage
        SNR = Total_TKE / Noise_TKE;
        noise_percentage = 100 * (1 / SNR);
        real_signal = 100 - noise_percentage;
        
        disp(['SNR: ', num2str(SNR), ', Signal: ', num2str(real_signal), '%']);
        
        % Calculate number of batches
        numBatches = ceil(setup.imProperties.imageCount / setup.imProperties.caseImages);
        
        % Calculate average cumulative energy across batches (one at a time)

        % --- Use domain for POD directory ---
        pod_dir = fullfile(stats_dir, domain, num2str(1));
         % Create output directory for reconstructed data
        output_dir = fullfile(data_dir, 'POD_Reconstructed');
        if ~exist(output_dir, 'dir')
            mkdir(output_dir);
        end
        
        % --- Save coordinates in output directory before processing images ---
        save(fullfile(output_dir, 'Co_ords.mat'), '-struct', 'Co_ords');
        [h, w] = size(y);
        Data_POD_current = load(fullfile(pod_dir, ['POD_stats_', num2str(w), 'x', num2str(h), '.mat']));
        
        cumulativeEnergy_avg = Data_POD_current.cumulativeEnergy;
        
        % Find required number of modes
        k_required = find(cumulativeEnergy_avg > real_signal, 1, 'first');
        if isempty(k_required)
            k_required = length(cumulativeEnergy_avg);
        end
        
        disp(['Modes required for ', num2str(real_signal), '% energy: ', num2str(k_required)]);
        
        % Get grid dimensions
        [h, w] = size(y);
        
       

        % Process each batch separately
        for batchNum = 1:numBatches
            disp(['Processing batch ', num2str(batchNum), ' of ', num2str(numBatches)]);
            
            % Load POD data for this batch only
            pod_dir = fullfile(stats_dir, domain, num2str(batchNum));
            if batchNum > 1
                Data_POD_current = load(fullfile(pod_dir, ['POD_stats_', num2str(w), 'x', num2str(h), '.mat']));
            end
            % Calculate image range for this batch
            start_img = (batchNum - 1) * setup.imProperties.caseImages + 1;
            end_img = min(batchNum * setup.imProperties.caseImages, setup.imProperties.imageCount);
            
            % Process all images in this batch
            parfor imNo = start_img:end_img
                kdx = imNo - (batchNum - 1) * setup.imProperties.caseImages;
                
                % Load original velocity data
                VelData = load(fullfile(data_dir, [num2str(sprintf(setup.instantaneous.nameConvention{1}, imNo))]));
                u_orig = VelData.piv_result(i).ux;
                v_orig = VelData.piv_result(i).uy;
                
                % Initialize reconstructed fields
                u_reconstructed = u_orig;
                v_reconstructed = v_orig;
                
                % --- Apply reconstruction based on domain ---
                UV_rec = Data_POD_current.U_svd(:, 1:k_required) * ...
                        (Data_POD_current.S_svd(1:k_required, 1:k_required) * ...
                         Data_POD_current.V_svd(kdx, 1:k_required)');
                
                u_rec = UV_rec(1:h*w) + Data_POD_current.mean_u;
                v_rec = UV_rec(h*w+1:end) + Data_POD_current.mean_v;
                
                u_rec = reshape(u_rec, [h, w]);
                v_rec = reshape(v_rec, [h, w]);
                
                if strcmp(domain, 'Full')
                    reconstruction_mask = ~mask;
                    u_reconstructed(reconstruction_mask) = u_rec(reconstruction_mask);
                    v_reconstructed(reconstruction_mask) = v_rec(reconstruction_mask);
                elseif strcmp(domain, 'Below')
                    u_reconstructed(cavity_region) = u_rec(cavity_region);
                    v_reconstructed(cavity_region) = v_rec(cavity_region);
                elseif strcmp(domain, 'Above')
                    above_region = (y > 0) & ~mask;
                    u_reconstructed(above_region) = u_rec(above_region);
                    v_reconstructed(above_region) = v_rec(above_region);
                end
                % ----------------------------------------------------------
                
                % Save reconstructed data
                piv_result = VelData.piv_result;
                piv_result(i).ux = u_reconstructed;
                piv_result(i).uy = v_reconstructed;
                piv_result(i).mask = mask; % Store the mask

                save_filename = fullfile(output_dir, [num2str(sprintf('%05d.mat', imNo))]);
                PIVSAVE(save_filename, piv_result);
                
                if mod(imNo, 1000) == 0
                    disp(['Processed image ', num2str(imNo), ' of ', num2str(setup.imProperties.imageCount)]);
                end
            end
            
            % Clear current batch data before loading next batch
            clear Data_POD_current;
            disp(['Completed batch ', num2str(batchNum)]);
        end
        
        % Save reconstruction parameters
        POD_params = struct();
        POD_params.modes_required = k_required;
        POD_params.SNR = SNR;
        POD_params.signal_percentage = real_signal;
        POD_params.noise_percentage = noise_percentage;
        POD_params.domain = domain;
        POD_params.uncertainty = uncertainty;
        POD_params.window_size = [setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)];
        
        save(fullfile(output_dir, 'POD_reconstruction_params.mat'), 'POD_params');
        
        disp(['POD reconstruction completed for run ', num2str(i)]);
        disp(['Output saved to: ', output_dir]);
    end
    end
end
