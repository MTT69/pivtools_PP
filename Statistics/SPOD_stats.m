function SPOD_stats(setup,run,CameraNo,endpoint,nDFT,combine)
    if setup.pipeline.SPOD
        directory_path = fullfile(setup.directory.base, 'CalibratedPIV', num2str(setup.imProperties.imageCount), ['Cam', num2str(CameraNo)], 'Instantaneous',endpoint);
        file_name = num2str(sprintf(setup.instantaneous.nameConvention{1}, 1));
        file_path = fullfile(directory_path, file_name);
        calibrated = load(file_path);
        
        [rows, cols] = size(calibrated.piv_result(run).ux); % Grid size
        U = zeros(setup.imProperties.imageCount,rows,cols ); % Preallocate u
        V = zeros(setup.imProperties.imageCount,rows,cols); % Preallocate v
        fprintf('Loading SPOD images at %s\n',(datetime('now')));
        parfor imNo = 1:setup.imProperties.imageCount
            % Construct the filename
            file_name = num2str(sprintf(setup.instantaneous.nameConvention{1}, imNo));
            file_path = fullfile(directory_path, file_name);
            calibrated = load(file_path);
            u= calibrated.piv_result(run).ux;
            v = calibrated.piv_result(run).uy;
            U(imNo,:,:) = u; % Store u component
            V(imNo,:,:) = v; % Store v component
            
        end
        if combine

            snapshot = [U,V];
            
        else
            snapshot = U;
        end
        clear U V
        nOvlp = nDFT / 2;
        temporalMean = (squeeze(mean(snapshot, 1)));
        fprintf('Computing SPOD statistics at %s\n',(datetime('now')));
        [L, P, f, ~, A] = spod(double(snapshot), nDFT, [], nOvlp, setup.imProperties.dt);

        directory = fullfile(directory_path,'SPOD');
        if ~exist(directory, 'dir')
            mkdir(directory);
        end
        % Define the file name and path for saving
        save_file_path = fullfile(directory, 'SPOD_results.mat');

        % Save all variables to the .mat file with version 7.3 format
        save(save_file_path, 'temporalMean', 'L', 'P', 'f', 'A', '-v7.3');

        % Assuming L is your energy array with size n x m (n frequencies, m modes)
        n = size(L, 1);  % Number of frequencies
        m = size(L, 2);  % Number of modes
        
        % Initialize a figure
        figure;
        fig = gcf;
        
        hold on;
        
        % Loop through each frequency (row) and compute cumulative sum for each
        for i = 1:n
            % Compute the cumulative sum for each frequency (row)
            cumulative_energy = cumsum(L(i, :));
            
            % Plot the cumulative energy for this frequency
            % Use the corresponding frequency value from f in the DisplayName
            plot(1:m, cumulative_energy, 'DisplayName', sprintf('f = %.2f', f(i)));
        end
        set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
        % Add labels and title
        xlabel('Mode');
        ylabel('Cumulative Modal Energy');
        title('Cumulative Modal Energy for Each Frequency');
        
        % Add legend
        legend show;
        saveas(fig, fullfile(directory, join(['Energy Breakdown' '.jpg'],"")));
        saveas(fig, fullfile(directory, join(['Energy Breakdown' '.epsc'],"")));
        saveas(fig, fullfile(directory, join(['Energy Breakdown' '.fig'],"")));
        close all
    end
    
   
    
  
   
   
    
    
    
    
end