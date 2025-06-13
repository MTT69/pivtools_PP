function RPCA_infil_top_bottom(setup, lambda, maxIter,tol, peakHeight, gappy, CameraNo, split)


    % RPCA_infil_top_bottom(setup, lambda, maxIter,tol, peakHeight, run, gappy, CameraNo, split)
        %
        % Performs Proper Orthogonal Decomposition (POD) on Particle Image Velocimetry 
        % (PIV) data to extract low-rank and sparse components using robust principal 
        % component analysis (RPCA). The function processes calibrated and uncalibrated
        % velocity field data, applying masks to exclude invalid vectors based on
        % specified criteria, and saves the results.
        %
        % Inputs:
        %   setup      - Struct containing the experimental setup and file structure:
        %                - setup.imProperties.imageCount: Number of PIV image pairs.
        %                - setup.instantaneous.nameConvention: File naming convention 
        %                  for instantaneous PIV results.
        %                - setup.directory.base: Base directory for file saving/loading.
        %   lambda     - Regularization parameter for RPCA. If set to 1, it defaults
        %                to 1/sqrt(imageCount).
        %   maxIter    - Maximum number of iterations for the RPCA algorithm.
        %   tol        - Convergence tolerance for RPCA.
        %   peakHeight - Threshold for peak magnitude in the uncalibrated PIV results.
        %   gappy      - A boolean that describes method of rpca reconstruction
        %   CameraNo   - Camera no identifier int
        %   split      - integer that is the row at which the top and bottom of the image are split. use RPCA_infil for full image
        %
        % Outputs:
        %   None (Results are saved directly to the file system).
        %
        % Description:
        %   - The function computes the low-rank (background) and sparse (fluctuation) 
        %     components of the PIV velocity fields (ux, uy) using RPCA.
        %   - Invalid velocity vectors are identified using masks based on 
        %     uncalibrated PIV results (e.g., peak choice and magnitude).
        %   - these vectors are set to zero for no energy contribution
        %   - The U and V vectors are ran as one large array - they are
        %     amplified by 10 x to ensure zero velocity mask stands out.
        %   - The gappy boolean means that only masked vectors are updated.
        %   - The results (low-rank and sparse components, rejected vectors) are 
        %     stored in the same structure as the original calibrated data and 
        %     saved to disk.
        %   - Parallel processing is utilized to handle large datasets efficiently.
        %
        % Example:
        %   % Perform POD decomposition on PIV data with default lambda for run 1
        %   POD_infil(setup, 1, 1000, 1e-6, 0.5, 1,1);
        %
        % Notes:
        %   - The RPCA algorithm is implemented using the `inexact_alm_rpc` function.
        %   - The function operates on a grid with size defined by the first calibrated 
        %     PIV result in the dataset.
        %   - Results are saved in the directory structure:
        %     'CalibratedPIV/<imageCount>/Cam<CameraNo>/Instantaneous/RPCA'.
        %
    
    
    if lambda ==0
        lambda = 1/sqrt(setup.imProperties.imageCount);
    end
    for i = 1:length(setup.instantaneous.runs)
        run = setup.instantaneous.runs(i);
    
        if setup.pipeline.RPCA_split
            directory_path = fullfile(setup.directory.base, 'CalibratedPIV', num2str(setup.imProperties.imageCount), ['Cam', num2str(CameraNo)], 'Instantaneous');
        
            file_name = num2str(sprintf(setup.instantaneous.nameConvention{1}, 1));
            file_path = fullfile(directory_path, file_name);
            calibrated = load(file_path);
            load(fullfile(directory_path, "Co_ords" ))
        
            [rows, cols] = size(calibrated.piv_result(run).ux); % Grid size
            U_top = zeros(split*cols, setup.imProperties.imageCount); % Preallocate u_top
            U_bottom = zeros((rows-split)*cols, setup.imProperties.imageCount); % Preallocate u_bottom
            V_top = zeros(split*cols, setup.imProperties.imageCount); % Preallocate v_top
            V_bottom = zeros((rows-split)*cols, setup.imProperties.imageCount); % Preallocate v_bottom
            
            masks_top = zeros(2*split*cols, setup.imProperties.imageCount);
            masks_bottom = zeros(2*(rows-split)*cols, setup.imProperties.imageCount);

            fprintf('Loading rpca images top_bottom at %s\n',(datetime('now')));
            Setup_parpool(setup, 'Images')
            for imNo = 1:setup.imProperties.imageCount
                % Construct the filename
                file_name = num2str(sprintf(setup.instantaneous.nameConvention{1}, imNo));
                file_path = fullfile(directory_path, file_name);
                calibrated = load(file_path);
                unCalibrated = load(fullfile(setup.directory.base, 'UncalibratedPIV', num2str(setup.imProperties.imageCount), ['Cam', num2str(1)], 'Instantaneous', sprintf( setup.instantaneous.nameConvention{1},imNo)));
                % Store u and v in respective arrays (assuming rectangular grids)
                mask = (unCalibrated.piv_result(run).peak_choice ~=1) | (unCalibrated.piv_result(run).peak_mag < peakHeight) | unCalibrated.piv_result(run).nan_mask ~=0;
                u= calibrated.piv_result(run).ux;
                u(mask) =0;
                v = calibrated.piv_result(run).uy;
                v(mask) =0;

                U_top(:,imNo) = reshape(u(1:split,:), [], 1); % Store u component
                U_bottom(:,imNo) = reshape(u(split+1:end,:), [], 1); % Store u component
                V_top(:,imNo) = reshape(v(1:split,:), [], 1); % Store v component
                V_bottom(:,imNo) = reshape(v(split+1:end,:), [], 1); % Store v component
                masks_top(:,imNo) = [reshape(mask(1:split,:), [], 1);reshape(mask(1:split,:), [], 1)];
                masks_bottom(:,imNo) = [reshape(mask(split+1:end,:), [], 1);reshape(mask(split+1:end,:), [], 1)];
            end
        
            % Reshape u and v into 2D matrices for decomposition
            % Rows: spatial points (setup.imProperties.imageCount*m); Columns: time instances
        
            U_top = U_top * 10;
            U_bottom = U_bottom * 10;
            V_top = V_top * 10;
            V_bottom = V_bottom * 10;
                    
            
            % Combine U and scaled V into X for top and bottom separately
            X_top = [U_top; V_top];
            X_bottom = [U_bottom; V_bottom];
            fprintf('Performing rpca top_bottom at %s\n',(datetime('now')));
            
            % Perform low-rank and sparse decomposition for top and bottom separately
            [L_top, S_top, ~] = inexact_alm_rpc(X_top, lambda, tol, maxIter, logical(masks_top), gappy);
            [L_bottom, S_bottom, ~] = inexact_alm_rpc(X_bottom, lambda, tol, maxIter, logical(masks_bottom), gappy);
            
            % Reshape decomposed matrices back to original dimensions for top and bottom separately
            Lu_top = reshape(L_top(1:end/2, :), split, cols, setup.imProperties.imageCount) / 10; % Low-rank u for top
            Lv_top = reshape(L_top(end/2+1:end, :), split, cols, setup.imProperties.imageCount) / 10; % Low-rank v for top
            Su_top = reshape(S_top(1:end/2, :), split, cols, setup.imProperties.imageCount) / 10; % Sparse u for top
            Sv_top = reshape(S_top(end/2+1:end, :), split, cols, setup.imProperties.imageCount) / 10; % Sparse v for top
            
            Lu_bottom = reshape(L_bottom(1:end/2, :), rows-split, cols, setup.imProperties.imageCount) / 10; % Low-rank u for bottom
            Lv_bottom = reshape(L_bottom(end/2+1:end, :), rows-split, cols, setup.imProperties.imageCount) / 10; % Low-rank v for bottom
            Su_bottom = reshape(S_bottom(1:end/2, :), rows-split, cols, setup.imProperties.imageCount) / 10; % Sparse u for bottom
            Sv_bottom = reshape(S_bottom(end/2+1:end, :), rows-split, cols, setup.imProperties.imageCount) / 10; % Sparse v for bottom
            rejected_mask_top = reshape(masks_top(1:end/2, :), split, cols, setup.imProperties.imageCount); % Rejected vectors for top
            rejected_mask_bottom = reshape(masks_bottom(1:end/2, :), rows-split, cols, setup.imProperties.imageCount); % Rejected vectors for bottom

            % Combine the top and bottom rejected masks
            rejected_mask = [rejected_mask_top; rejected_mask_bottom];

            fprintf('Saving rpca vectors at %s\n',(datetime('now')));
            directory = fullfile(directory_path,'RPCA_top_bottom');
            if ~exist(directory, 'dir')
                mkdir(directory);
            end
            save(fullfile(directory,"Co_ords"),"Co_ords")
            Setup_parpool(setup, 'Images')
            for imNo = 1:setup.imProperties.imageCount
                % Construct the filename
                
                file_name = num2str(sprintf(setup.instantaneous.nameConvention{1}, imNo));
                file_path = fullfile(directory, file_name);
                unCalibrated = load(fullfile(setup.directory.base, 'UncalibratedPIV', num2str(setup.imProperties.imageCount), ['Cam', num2str(1)], 'Instantaneous', sprintf( setup.instantaneous.nameConvention{1},imNo)));
                b_mask = unCalibrated.piv_result(run).b_mask;
        
                piv_result=struct();
                % Rebuild the full ux, uy, and su by combining the top and bottom
                ux_full = [Lu_top(:,:,imNo); Lu_bottom(:,:,imNo)];
                uy_full = [Lv_top(:,:,imNo); Lv_bottom(:,:,imNo)];
                su_full = [Su_top(:,:,imNo); Su_bottom(:,:,imNo)];
                sv_full = [Sv_top(:,:,imNo); Sv_bottom(:,:,imNo)];

                piv_result(run).b_mask = b_mask;
                piv_result(run).ux = ux_full;
                piv_result(run).uy = uy_full;
                piv_result(run).ux_s = su_full;
                piv_result(run).uy_s = sv_full;
                piv_result(run).Rejected_vectors = rejected_mask(:,:,imNo);

                PIVSAVE(file_path,piv_result)
        
                
            end

            ux_full = [Lu_top; Lu_bottom];
            uy_full = [Lv_top; Lv_bottom];

            % Calculate means of ux and uy
            mean_ux = reshape(mean(ux_full, 3), rows, cols);
            mean_uy = reshape(mean(uy_full, 3), rows, cols);

            % Save means to file
            % Define the file path
            mean_stats_run = fullfile(setup.directory.base, 'Statistics', ...
                num2str(setup.imProperties.imageCount), ...
                ['Cam' num2str(CameraNo)], 'Instantaneous', ...
                'Calibrated', 'RPCA_top_bottom', ...
                ['mean_stats_run_', num2str(run), '.mat']);
            
            % Extract the folder path
            folder_path = fileparts(mean_stats_run);
            
            % Ensure the folder exists
            if ~isfolder(folder_path)
                mkdir(folder_path);
            end


            save(mean_stats_run, 'mean_ux', 'mean_uy');
            
        
            fprintf('finished at %s\n',(datetime('now')));
        end
    end
end