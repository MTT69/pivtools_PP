function plot_SPOD_reconstruction(setup, CameraNo, endpoint, combine, modes, frequencies, numTimeElements,run)
    if setup.pipeline.SPOD
        % Dynamically construct the save directory path
        directory = fullfile(setup.directory.base, 'CalibratedPIV', num2str(setup.imProperties.imageCount), ...
            ['Cam', num2str(CameraNo)], 'Instantaneous', endpoint, 'SPOD');
        
        SPOD_stats = load(fullfile(directory, 'SPOD_results.mat'));
        save_directory = fullfile(directory, ...
            ['Reconstruction modes ' num2str(modes(1)), ' _ ' num2str(modes(end)), ' Frequencies ' num2str(SPOD_stats.f(frequencies(1))), ' _ ', num2str(SPOD_stats.f(frequencies(2))) 'hz']);
        
        % Ensure the directory exists
        if ~exist(save_directory, 'dir')
            error('Save directory does not exist: %s', save_directory);
        end
        
        % Prepare to create the GIF
        gif_filename = fullfile(save_directory, 'FlowDevelopment.gif');
        first_frame = true;  % Flag to handle the first frame for the GIF
        
        % Loop through the specified number of time elements
        for imNo = 1:numTimeElements
            % Load the current PIV result file
            file_name = num2str(sprintf(setup.instantaneous.nameConvention{1}, imNo));
            filePath = fullfile(save_directory, file_name);
            
            data = load(filePath);
            
            
            % Extract the PIV result
            piv_result = data.piv_result(run);
            
            % Create a new figure
            fig = figure;
            
            if combine
                % Handle combined mode
                if isfield(piv_result, 'ux') && isfield(piv_result, 'uv')
                    % Plot the ux (x-velocity)
                    subplot(2, 1, 1);
                    imagesc(piv_result.ux);
                    daspect([1 1 1]);
                    colorbar;
                    title(sprintf('U Velocity - Time Frame %d', imNo));
                    xlabel('X Coordinate');
                    ylabel('Y Coordinate');
                    
                    % Plot the uv (y-velocity)
                    subplot(2, 1, 2);
                    imagesc(piv_result.uv);
                    daspect([1 1 1]);
                    colorbar;
                    title(sprintf('V Velocity - Time Frame %d', imNo));
                    xlabel('X Coordinate');
                    ylabel('Y Coordinate');
                else
                    warning('File %s does not contain both "ux" and "uv". Skipping...', file_name{imNo});
                    close(fig);
                    continue;
                end
            else
                % Handle non-combine mode
                if isfield(piv_result, 'ux')
                    % Create the plot for 'ux' component
                    imagesc(piv_result.ux);
                    daspect([1 1 1]);
                    colorbar;
                    title(sprintf('Flow Development - U Velocity - Time Frame %d', imNo));
                    xlabel('X Coordinate');
                    ylabel('Y Coordinate');
                else
                    warning('File %s does not contain "ux". Skipping...', file_name{imNo});
                    close(fig);
                    continue;
                end
            end
            
            % Capture the current frame
            frame = getframe(fig);
            
            % Capture the current frame

            
            % If it's the first frame, initialize the GIF
            if first_frame
                % Convert the frame to the correct format for saving as GIF
                [A, map] = rgb2ind(frame.cdata, 256);  % Convert to indexed image with a colormap
                
                % Save the first frame as the GIF
                imwrite(A, map, gif_filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.1);
                first_frame = false;  % No longer the first frame
            else
                % Convert the frame to the correct format for saving as GIF
                [A, map] = rgb2ind(frame.cdata, 256);  % Convert to indexed image with a colormap
                
                % Append the frame to the GIF
                imwrite(A, map, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
            end
            
            % Close the current figure after saving the frame
            close(fig);

        
        end
    end
end
