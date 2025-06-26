function Gamma2_video(setup, endpoint, CameraNo, type, use_merged, acq_freq)
% gamma2_VIDEO - Function to create gamma2 field visualization video
%   This function generates a video showing instantaneous gamma2 field
%
%   Usage: Gamma2_video(setup, endpoint, CameraNo, type, use_merged, acq_freq)
%   type: 'Instantaneous', 'Calibrated', or 'Merged'
%   acq_freq: acquisition frequency in Hz

if nargin < 4
    type = 'Instantaneous'; % Default for backward compatibility
end
if nargin < 5
    use_merged = false;
end
if nargin < 6
    error('Acquisition frequency (acq_freq) must be provided.');
end

if setup.pipeline.gamma2_video
    fprintf('Creating gamma2 video for base: %s at %s\n', setup.directory.base, datetime('now'));
    
    % Get data paths using centralized function
    paths = get_data_paths(setup, type, endpoint, CameraNo, use_merged);
    dataloc = paths.data_dir;
    output_dir = fullfile(paths.video_dir, 'gamma2_Analysis');
    
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    for i = setup.instantaneous.runs
        fprintf('Processing run %d for gamma2 video...\n', i);
        
        % Load coordinates
        Co_ords = load(fullfile(dataloc, 'Co_ords.mat'));
        x = Co_ords.Co_ords(i).x;
        y = Co_ords.Co_ords(i).y;
        
        % Get spatial extent
        xcorners = [x(1,1), x(end, end)];
        ycorners = [y(1,1), y(end, end)];
        
        % Load first image to get mask
        VelData = load(fullfile(dataloc, sprintf('%05d.mat', 1)));
        b_mask = VelData.piv_result(i).b_mask;
        
        % Pre-calculate gamma2 limits for consistent colorbar scaling
        fprintf('Calculating gamma2 limits for consistent scaling...\n');
        sample_frames = 1:50;
        all_gamma2_values = [];
        
        for sample_idx = 1:length(sample_frames)
            sample_imNo = sample_frames(sample_idx);
            VelData_sample = load(fullfile(dataloc, sprintf('%05d.mat', sample_imNo)));
            u_sample = VelData_sample.piv_result(i).ux;
            v_sample = VelData_sample.piv_result(i).uy;
            
            u_smooth = u_sample;
            v_smooth = v_sample;
            
            % Apply mask to smoothed data
            u_smooth(b_mask) = NaN;
            v_smooth(b_mask) = NaN;
            
            % Calculate gamma2 field
            gammaField_sample = gamma2(x, y, u_smooth, v_smooth, 10);
            gammaField_sample(b_mask) = NaN;
            
            valid_gamma2 = gammaField_sample(~isnan(gammaField_sample));
            all_gamma2_values = [all_gamma2_values; valid_gamma2(:)];
            
            if mod(sample_idx, 10) == 0
                fprintf('Sampled %d/%d frames for limit calculation\n', sample_idx, length(sample_frames));
            end
        end
        
        % Calculate consistent limits using percentiles
        gamma2_lower_limit = prctile(all_gamma2_values, 5);
        gamma2_upper_limit = prctile(all_gamma2_values, 95);
        
        % Ensure symmetric limits around zero for better visualization
        max_abs_limit = max(abs(gamma2_lower_limit), abs(gamma2_upper_limit));
        gamma2_lower_limit = -max_abs_limit;
        gamma2_upper_limit = max_abs_limit;
        
        fprintf('gamma2 limits: [%.2f, %.2f]\n', gamma2_lower_limit, gamma2_upper_limit);
        
        % Video setup
        video_filename = fullfile(output_dir, sprintf('gamma2_Video_Run%d_%dx%d.mp4', ...
                                 i, setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)));
        v = VideoWriter(video_filename, 'MPEG-4');
        v.FrameRate = 30; % Adjust as needed
        open(v);
        
        % Process each frame
        for imNo = 1:setup.imProperties.caseImages
            % Load velocity data
            VelData = load(fullfile(dataloc, sprintf('%05d.mat', imNo)));
            u_vel = VelData.piv_result(i).ux;
            v_vel = VelData.piv_result(i).uy;
            
            
            % Apply mask to smoothed data
            u_vel(b_mask) = NaN;
            v_vel(b_mask) = NaN;
            
            % Calculate gamma2 field using smoothed data
            gammaField = gamma2(x, y, u_vel, v_vel, 10);
            gammaField(b_mask) = 0;
            
            % Create figure using plot_save_mask with return_figure option
            variableName = sprintf('gamma2 Field');
            
            % Get figure handle instead of saving files
            fig = plot_save_mask(gammaField, b_mask, xcorners, ycorners, ...
                          setup.figures.axisFontSize, setup.figures.titleFontSize, ...
                          output_dir, setup.instantaneous.runs, ...
                          [setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)], ...
                          i, variableName, 'Inst', output_dir, ...
                          gamma2_lower_limit, gamma2_upper_limit, true);
            
            set(fig, 'Visible', 'off');
            set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

            
            % Get current axes
            ax = gca;
            
            % Add xlabel and ylabel
            xlabel(ax, 'x [mm]', 'Interpreter', 'latex', 'FontSize', setup.figures.labelFontSize);
            ylabel(ax, 'y [mm]', 'Interpreter', 'latex', 'FontSize', setup.figures.labelFontSize);
            
            % Update colorbar label
            cb = colorbar(ax);
            cb.Label.String = 'Gamma2';
            cb.Label.Interpreter = 'latex';
            cb.Label.FontSize = setup.figures.labelFontSize;
            
            % Add time to title
            time_s = (imNo - 1) / acq_freq;
            title(ax, sprintf('gamma2 Field, Time: %.3fs', time_s), ...
                  'Interpreter', 'latex', 'FontSize', setup.figures.titleFontSize);
            
            % Capture frame
            frame = getframe(fig);
            writeVideo(v, frame);
            
            % Close figure
            close(fig);
            
            % Progress indicator
            if mod(imNo, 10) == 0
                fprintf('Processed frame %d/%d\n', imNo, setup.imProperties.caseImages);
            end
        end
        
        % Close video
        close(v);
        fprintf('Video saved: %s\n', video_filename);
        
       
    end
    
    fprintf('gamma2 video generation completed at %s\n', datetime('now'));
end
