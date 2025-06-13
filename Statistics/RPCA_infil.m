function RPCA_infil(setup, lambda, maxIter,tol, peakHeight,gappy,CameraNo)


% RPCA_infil(setup, lambda, maxIter,tol, peakHeight,run,gappy,CameraNo)
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
    %   run        - Index of the current PIV run (used for loading/saving results).
    %   gappy      - A boolean that describes method of rpca reconstruction
    %   CameraNo   - Camera no identifier int
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
    %     amplified by 10 as a baseline and then U and V are normalised to 
    %     ensure zero velocity mask stands out and that the error of U and V are considered the same.
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

    if setup.pipeline.RPCA
        directory_path = fullfile(setup.directory.base, 'CalibratedPIV', num2str(setup.imProperties.imageCount), ['Cam', num2str(CameraNo)], 'Instantaneous');

        file_name = num2str(sprintf(setup.instantaneous.nameConvention{1}, 1));
        file_path = fullfile(directory_path, file_name);
        calibrated = load(file_path);
        load(fullfile(directory_path, "Co_ords" ))

        [rows, cols] = size(calibrated.piv_result(run).ux); % Grid size
        U = zeros(rows* cols, setup.imProperties.imageCount); % Preallocate u
        V = zeros(rows* cols, setup.imProperties.imageCount); % Preallocate v
        masks = zeros(2*rows* cols, setup.imProperties.imageCount);
        fprintf('Loading rpca images at %s\n',(datetime('now')));
        for imNo = 1:setup.imProperties.imageCount
            % Construct the filename
            file_name = num2str(sprintf(setup.instantaneous.nameConvention{1}, imNo));
            file_path = fullfile(directory_path, file_name);
            calibrated = load(file_path);
            unCalibrated = load(fullfile(setup.directory.base, 'UncalibratedPIV', num2str(setup.imProperties.imageCount), ['Cam', num2str(1)], 'Instantaneous', sprintf( setup.instantaneous.nameConvention{1},imNo)));
            % Store u and v in respective arrays (assuming rectangular grids)
            mask = (unCalibrated.piv_result(run).peak_choice ~=1) | (unCalibrated.piv_result(run).peak_mag < peakHeight);
            u= calibrated.piv_result(run).ux;
            u(mask) =0;
            v = calibrated.piv_result(run).uy;
            v(mask) =0;
            U(:,imNo) = u(:); % Store u component
            V(:,imNo) = v(:); % Store v component
            masks(:,imNo) = [mask(:);mask(:)];
        end

        % Reshape u and v into 2D matrices for decomposition
        % Rows: spatial points (setup.imProperties.imageCount*m); Columns: time instances
        U = U*10;
        V = V*10;
        
        % Combine U and scaled V into X
        X = [U; V];
        fprintf('Performing rpca at %s\n',(datetime('now')));

        % Perform low-rank and sparse decomposition
        [L, S, ~] = inexact_alm_rpc(X, lambda, tol, maxIter,logical(masks),gappy);

        % Reshape decomposed matrices back to original dimensions
        Lu = reshape(L(1:end/2, :), rows, cols, setup.imProperties.imageCount)/10; % Low-rank u
        Lv = reshape(L(end/2+1:end, :), rows, cols, setup.imProperties.imageCount)/10; % Low-rank v
        Su = reshape(S(1:end/2, :), rows, cols, setup.imProperties.imageCount)/10; % Sparse u
        Sv = reshape(S(end/2+1:end, :), rows, cols, setup.imProperties.imageCount)/10; % Sparse v
        rejected_mask = reshape(masks(1:end/2,:),rows,cols,setup.imProperties.imageCount);
        fprintf('Saving rpca images at %s\n',(datetime('now')));
        directory = fullfile(directory_path,'RPCA');
        if ~exist(directory, 'dir')
            mkdir(directory);
        end
        save(fullfile(directory,"Co_ords"),"Co_ords")
        for imNo = 1:setup.imProperties.imageCount
            % Construct the filename
            
            file_name = num2str(sprintf(setup.instantaneous.nameConvention{1}, imNo));
            file_path = fullfile(directory, file_name);
            unCalibrated = load(fullfile(setup.directory.base, 'UncalibratedPIV', num2str(setup.imProperties.imageCount), ['Cam', num2str(1)], 'Instantaneous', sprintf( setup.instantaneous.nameConvention{1},imNo)));
            b_mask = unCalibrated.piv_result(run).b_mask;

            piv_result=struct();
            piv_result(run).b_mask = b_mask;
            piv_result(run).ux = Lu(:,:,imNo);
            piv_result(run).uy = Lv(:,:,imNo);
            piv_result(run).ux_s = Su(:,:,imNo);
            piv_result(run).uy_s = Sv(:,:,imNo);
            piv_result(run).Rejected_vectors = rejected_mask(:,:,imNo);
            PIVSAVE(file_path,piv_result)

            
        end
        mean_ux = reshape(mean(Lu, 3), rows, cols);
        mean_uy = reshape(mean(Lv, 3), rows, cols);
        mean_stats_run = fullfile(setup.directory.base, 'Statistics', num2str(setup.imProperties.imageCount), ['Cam' num2str(CameraNo)], 'Instantaneous', 'Calibrated', 'RPCA', ['mean_stats_run_', num2str(run), '.mat']);
        save(mean_stats_run, 'mean_ux', 'mean_uy');

        fprintf('finished at %s\n',(datetime('now')));
    end
end
end