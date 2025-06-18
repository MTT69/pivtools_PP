function Pressure_video(setup, endpoint, CameraNo, type, use_merged)
% PRESSURE_VIDEO - Function to create pressure field visualization video
%   Usage: Pressure_video(setup, endpoint, CameraNo, type, use_merged)

if nargin < 4
    type = 'Instantaneous'; % Default to Instantaneous for backward compatibility
end
if nargin < 5
    use_merged = false; % Default to traditional camera data
end

if setup.pipeline.pressure_video
    fprintf('Creating pressure video for base: %s at %s\n', setup.directory.base, datetime('now'));
    
    % Get data paths using centralized function
    paths = get_data_paths(setup, type, endpoint, CameraNo, use_merged);
    dataloc = paths.data_dir;
    pressure_dataloc = paths.pressure_dir;
    output_dir = fullfile(paths.video_dir, 'Pressure_Analysis');
    
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % Check if pressure data exists
    if ~exist(pressure_dataloc, 'dir')
        error('Pressure data not found. Please run pressure_reconstruction first.');
    end
    
    for i = setup.instantaneous.runs
        fprintf('Processing run %d for pressure video...\n', i);
        
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
        
        % Pre-calculate pressure limits for consistent colorbar scaling
        fprintf('Calculating pressure limits for consistent scaling...\n');
        sample_frames = 1:50:setup.imProperties.imageCount; % Sample every 50th frame
        all_pressure_values = [];
        
        for sample_idx = 1:length(sample_frames)
            sample_imNo = sample_frames(sample_idx);
            try
                PressureData = load(fullfile(pressure_dataloc, sprintf('%05d.mat', sample_imNo)));
                pressure_sample = PressureData.pressure_result.P;
                
                % Apply mask
                pressure_sample(b_mask) = NaN;
                
                valid_pressure = pressure_sample(~isnan(pressure_sample));
                all_pressure_values = [all_pressure_values; valid_pressure(:)];
                
                if mod(sample_idx, 10) == 0
                    fprintf('Sampled %d/%d frames for limit calculation\n', sample_idx, length(sample_frames));
                end
            catch
                fprintf('Warning: Could not load pressure data for frame %d\n', sample_imNo);
                continue;
            end
        end
        
        if isempty(all_pressure_values)
            error('No valid pressure data found for limit calculation');
        end
        
        % Calculate consistent limits using percentiles
        pressure_lower_limit = prctile(all_pressure_values, 5);
        pressure_upper_limit = prctile(all_pressure_values, 95);
        
        % Ensure symmetric limits around zero for better visualization
        max_abs_limit = max(abs(pressure_lower_limit), abs(pressure_upper_limit));
%         pressure_lower_limit = -max_abs_limit;
%         pressure_upper_limit = max_abs_limit;
        pressure_lower_limit =-5;
        pressure_upper_limit =5;
        
        fprintf('Pressure limits: [%.2f, %.2f] Pa\n', pressure_lower_limit, pressure_upper_limit);
        
        % Video setup
        video_filename = fullfile(output_dir, sprintf('Pressure_Video_Run%d_%dx%d.mp4', ...
                                 i, setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)));
        v = VideoWriter(video_filename, 'MPEG-4');
        v.FrameRate = 30; % Adjust as needed
        open(v);
        
        % Process each frame
        for imNo = 1:setup.imProperties.imageCount
            try
                % Load pressure data
                PressureData = load(fullfile(pressure_dataloc, sprintf('%05d.mat', imNo)));
                pressure_field = PressureData.pressure_result.P;
                
                % Create figure using plot_save_mask with return_figure option
                variableName = sprintf('Instantaneous Pressure Field');
                
                % Get figure handle instead of saving files
                fig = plot_save_mask(pressure_field, b_mask, xcorners, ycorners, ...
                              setup.figures.axisFontSize, setup.figures.titleFontSize, ...
                              output_dir, setup.instantaneous.runs, ...
                              [setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)], ...
                              i, variableName, 'Inst', output_dir, ...
                              pressure_lower_limit, pressure_upper_limit, true);
                
                set(fig, 'Visible', 'off');
                
                % Get current axes
                ax = gca;
                
                % Update title with frame information
                title(ax, 'Instantaneous Pressure Field', ...
                      'Interpreter', 'latex', 'FontSize', setup.figures.titleFontSize);
                
                % Add xlabel and ylabel
                xlabel(ax, 'x [mm]', 'Interpreter', 'latex', 'FontSize', setup.figures.labelFontSize);
                ylabel(ax, 'y [mm]', 'Interpreter', 'latex', 'FontSize', setup.figures.labelFontSize);
                
                % Update colorbar label
                cb = colorbar(ax);
                cb.Label.String = 'Pressure [Pa]';
                cb.Label.Interpreter = 'latex';
                cb.Label.FontSize = setup.figures.labelFontSize;
                
                % Capture frame
                frame = getframe(fig);
                writeVideo(v, frame);
                
                % Close figure
                close(fig);
                
            catch ME
                fprintf('Warning: Could not process frame %d: %s\n', imNo, ME.message);
                continue;
            end
            
            % Progress indicator
            if mod(imNo, 100) == 0
                fprintf('Processed frame %d/%d\n', imNo, setup.imProperties.imageCount);
            end
        end
        
        % Close video
        close(v);
        fprintf('Video saved: %s\n', video_filename);
        
        % Create a sample frame for reference
        sample_frame_dir = fullfile(output_dir, 'Sample_Frames');
        if ~exist(sample_frame_dir, 'dir')
            mkdir(sample_frame_dir);
        end
        
        % Generate sample frame from middle of sequence
        mid_frame = round(setup.imProperties.imageCount / 2);
        try
            PressureData = load(fullfile(pressure_dataloc, sprintf('%05d.mat', mid_frame)));
            pressure_field = PressureData.pressure_result.P;
            
            % Create sample frame using traditional save method
            sample_variableName = sprintf('SampleFrame_Pressure_Run%d', i);
            
            plot_save_mask(pressure_field, b_mask, xcorners, ycorners, ...
                          setup.figures.axisFontSize, setup.figures.titleFontSize, ...
                          sample_frame_dir, setup.instantaneous.runs, ...
                          [setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)], ...
                          i, sample_variableName, 'Inst', sample_frame_dir, ...
                          pressure_lower_limit, pressure_upper_limit);
            
            % Load and enhance sample frame
            sample_fig_filename = fullfile(sample_frame_dir, [sample_variableName, num2str(setup.instantaneous.windowSize(i,1)), 'x', num2str(setup.instantaneous.windowSize(i,2)), '.fig']);
            fig = openfig(sample_fig_filename);
            ax = gca;
            
            title(ax, sprintf('Sample Frame: Instantaneous Pressure Field\\nFrame %d, Run %d', ...
                         mid_frame, i), ...
                  'Interpreter', 'latex', 'FontSize', setup.figures.titleFontSize);
            
            xlabel(ax, 'x [mm]', 'Interpreter', 'latex', 'FontSize', setup.figures.labelFontSize);
            ylabel(ax, 'y [mm]', 'Interpreter', 'latex', 'FontSize', setup.figures.labelFontSize);
            
            cb = colorbar(ax);
            cb.Label.String = 'Pressure [Pa]';
            cb.Label.Interpreter = 'latex';
            cb.Label.FontSize = setup.figures.labelFontSize;
            
            % Save enhanced sample frame
            enhanced_base = sprintf('SampleFrame_Pressure_Run%d_%dx%d_enhanced', i, setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2));
            saveas(fig, fullfile(sample_frame_dir, [enhanced_base, '.fig']));
            saveas(fig, fullfile(sample_frame_dir, [enhanced_base, '.png']));
            saveas(fig, fullfile(sample_frame_dir, [enhanced_base, '.eps']));
            close(fig);
            
        catch ME
            fprintf('Warning: Could not create sample frame: %s\n', ME.message);
        end
    end
    
    fprintf('Pressure video generation completed at %s\n', datetime('now'));
end

