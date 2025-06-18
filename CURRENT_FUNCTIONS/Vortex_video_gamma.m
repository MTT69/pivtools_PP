function Vortex_video_gamma(setup, endpoint, CameraNo, type, use_merged)
% VORTEX_VIDEO_GAMMA - Function to create vortex visualization video
%   This function generates a video showing instantaneous vorticity with gamma1 contour overlays
%
%   Usage: Vortex_video_gamma(setup, endpoint, CameraNo, type, use_merged)
%   type: 'Instantaneous', 'Calibrated', or 'Merged'
%   use_merged: true for merged data, false for traditional camera data

if nargin < 4
    type = 'Instantaneous'; % Default for backward compatibility
end
if nargin < 5
    use_merged = false; % Default to traditional camera data
end

if setup.pipeline.vortex_video_gamma
    fprintf('Creating vortex video for base: %s at %s\n', setup.directory.base, datetime('now'));
    
    % Get data paths using centralized function
    paths = get_data_paths(setup, type, endpoint, CameraNo, use_merged);
    dataloc = paths.data_dir;
    output_dir = fullfile(paths.video_dir, 'Vortex_Analysis');
    
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    for i = setup.instantaneous.runs
        fprintf('Processing run %d for vortex video...\n', i);
        
        % Load coordinates
        Co_ords = load(fullfile(dataloc, 'Co_ords.mat'));
        x = Co_ords.Co_ords(i).x;
        y = Co_ords.Co_ords(i).y;
        
        % Calculate gradients for vorticity calculation
        dx = gradient(x);
        [~, dy] = gradient(y);
        dy = -dy; % Flip y gradient
        
        % Get spatial extent
        xcorners = [x(1,1), x(end, end)];
        ycorners = [y(1,1), y(end, end)];
        
        % Load first image to get mask
        VelData = load(fullfile(dataloc, sprintf('%05d.mat', 1)));
        b_mask = VelData.piv_result(i).b_mask;
        
        % Pre-calculate vorticity limits for consistent colorbar scaling
        fprintf('Calculating vorticity limits for consistent scaling...\n');
        sample_frames = 1:50:setup.imProperties.imageCount; % Sample every 50th frame
        all_vorticity_values = [];
        
        for sample_idx = 1:length(sample_frames)
            sample_imNo = sample_frames(sample_idx);
            VelData_sample = load(fullfile(dataloc, sprintf('%05d.mat', sample_imNo)));
            u_sample = VelData_sample.piv_result(i).ux;
            v_sample = VelData_sample.piv_result(i).uy;
            
            [dvx_sample, ~] = gradient(v_sample);
            [~, duy_sample] = gradient(u_sample);
            duy_sample = -duy_sample;
            vorticity_sample = dvx_sample./dx - duy_sample./dy;
            vorticity_sample(b_mask) = NaN;
            
            valid_vorticity = vorticity_sample(~isnan(vorticity_sample));
            all_vorticity_values = [all_vorticity_values; valid_vorticity(:)];
            
            if mod(sample_idx, 10) == 0
                fprintf('Sampled %d/%d frames for limit calculation\n', sample_idx, length(sample_frames));
            end
        end
        
        % Calculate consistent limits using percentiles
        vort_lower_limit = prctile(all_vorticity_values, 5);
        vort_upper_limit = prctile(all_vorticity_values, 95);
        
        % Ensure symmetric limits around zero for better visualization
        max_abs_limit = max(abs(vort_lower_limit), abs(vort_upper_limit));
        vort_lower_limit = -max_abs_limit;
        vort_upper_limit = max_abs_limit;
%         vort_lower_limit = -0.005;
%         vort_upper_limit = 0.005;
        
        fprintf('Vorticity limits: [%.2f, %.2f]\n', vort_lower_limit, vort_upper_limit);
        
        % Video setup
        video_filename = fullfile(output_dir, sprintf('Vortex_Video_Run%d_%dx%d.mp4', ...
                                 i, setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)));
        v = VideoWriter(video_filename, 'MPEG-4');
        v.FrameRate = 30; % Adjust as needed
        open(v);
        
        % Process each frame
        for imNo = 1:setup.imProperties.imageCount
            % Load velocity data
            VelData = load(fullfile(dataloc, sprintf('%05d.mat', imNo)));
            u = VelData.piv_result(i).ux;
            v_vel = VelData.piv_result(i).uy;
            
            % Apply spatial smoothing to reduce noise
            % Option 1: Gaussian smoothing
            sigma = 1; % Adjust smoothing strength (higher = more smoothing)
%             u_smooth = imgaussfilt(u, sigma, 'FilterDomain', 'spatial');
%             v_smooth = imgaussfilt(v_vel, sigma, 'FilterDomain', 'spatial');
            u_smooth = u;
            v_smooth = v;
            
            
            % Apply mask to smoothed data
            u_smooth(b_mask) = NaN;
            v_smooth(b_mask) = NaN;
            
            % Calculate instantaneous vorticity using smoothed data
            [dvx, ~] = gradient(v_smooth);
            [~, duy] = gradient(u_smooth);
            duy = -duy;
            vorticity = dvx./dx - duy./dy;
            
            % Calculate gamma1 field for vortex detection using smoothed data
            gammaField = gamma1(x, y, u_smooth, v_smooth, 10);
            gammaField(b_mask) = 0;
            
            % Create figure using plot_save_mask with return_figure option
            variableName = sprintf('Vorticity with Gamma Contours');
            
            % Get figure handle instead of saving files
            fig = plot_save_mask(vorticity, b_mask, xcorners, ycorners, ...
                          setup.figures.axisFontSize, setup.figures.titleFontSize, ...
                          output_dir, setup.instantaneous.runs, ...
                          [setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)], ...
                          i, variableName, 'Inst', output_dir, ...
                          vort_lower_limit, vort_upper_limit, true);
            
            set(fig, 'Visible', 'off');
            
            % Get current axes and add gamma1 contours
            ax = gca;
            hold(ax, 'on');
            
            % Add positive and negative gamma1 contours with different styles
            [C_pos, h_pos] = contour(ax, x, y, gammaField, [0.5, 1, 2], 'LineWidth', 2, 'LineStyle', '-');
            [C_neg, h_neg] = contour(ax, x, y, gammaField, [-2, -1, -0.5], 'LineWidth', 2, 'LineStyle', '--');
            h_pos.LineColor = 'yellow';
            h_neg.LineColor = 'cyan';
            
            % Add vector overlay with reduced density and smoothed vectors
            vector_skip = 4; % Show every 4th vector (reduce density)
            [X_vec, Y_vec] = meshgrid(1:vector_skip:size(x,2), 1:vector_skip:size(x,1));
            
            % Subsample coordinates and velocities
            x_sub = x(1:vector_skip:end, 1:vector_skip:end);
            y_sub = y(1:vector_skip:end, 1:vector_skip:end);
            u_sub = u_smooth(1:vector_skip:end, 1:vector_skip:end);
            v_sub = v_smooth(1:vector_skip:end, 1:vector_skip:end);
            mask_sub = b_mask(1:vector_skip:end, 1:vector_skip:end);
            
            % Remove vectors in masked regions
            u_sub(mask_sub) = NaN;
            v_sub(mask_sub) = NaN;
            
            % Scale vectors for better visibility
            vector_scale = 0.8; % Adjust scale factor
            quiver(ax, x_sub, y_sub, u_sub*vector_scale, v_sub*vector_scale, 0, ...
                   'Color', 'white', 'LineWidth', 1.2, 'MaxHeadSize', 0.3, ...
                   'AutoScale', 'off');
            
%             % Update title with frame information
%             title(ax, sprintf('Instantaneous Vorticity with $\\Gamma_1$ Contours and Vectors\\nFrame %d/%d, Run %d', ...
%                          imNo, setup.imProperties.imageCount, i), ...
%                   'Interpreter', 'latex', 'FontSize', setup.figures.titleFontSize);
            
            % Add xlabel and ylabel
            xlabel(ax, 'x [mm]', 'Interpreter', 'latex', 'FontSize', setup.figures.labelFontSize);
            ylabel(ax, 'y [mm]', 'Interpreter', 'latex', 'FontSize', setup.figures.labelFontSize);
            
            % Add legend for contours
            legend_entries = {
                plot(ax, NaN, NaN, 'y-', 'LineWidth', 2), 'Positive $\Gamma_1$', ...
                plot(ax, NaN, NaN, 'c--', 'LineWidth', 2), 'Negative $\Gamma_1$'
            };
            legend(ax, [legend_entries{1}, legend_entries{3}], {legend_entries{2}, legend_entries{4}}, ...
                   'Location', 'northeast', 'Interpreter', 'latex', ...
                   'FontSize', setup.figures.legendFontSize);
            
            % Update colorbar label
            cb = colorbar(ax);
            cb.Label.String = 'Vorticity [1/s]';
            cb.Label.Interpreter = 'latex';
            cb.Label.FontSize = setup.figures.labelFontSize;
            
            hold(ax, 'off');
            
            % Capture frame
            frame = getframe(fig);
            writeVideo(v, frame);
            
            % Close figure - no file cleanup needed since we didn't save files
            close(fig);
            
            % Progress indicator
            if mod(imNo, 10) == 0
                fprintf('Processed frame %d/%d\n', imNo, setup.imProperties.imageCount);
            end
        end
        
        % Close video
        close(v);
        fprintf('Video saved: %s\n', video_filename);
        
        % Create a sample frame for reference using the same approach
        sample_frame_dir = fullfile(output_dir, 'Sample_Frames');
        if ~exist(sample_frame_dir, 'dir')
            mkdir(sample_frame_dir);
        end
        
        mid_frame = round(setup.imProperties.imageCount / 2);
        VelData = load(fullfile(dataloc, sprintf('%05d.mat', mid_frame)));
        u = VelData.piv_result(i).ux;
        v_vel = VelData.piv_result(i).uy;
        
        [dvx, ~] = gradient(v_vel);
        [~, duy] = gradient(u);
        duy = -duy;
        vorticity = dvx./dx - duy./dy;
        gammaField = gamma1(x, y, u, v_vel, 10);
        gammaField(b_mask) = 0;
        
        % Create sample frame using traditional save method
        sample_variableName = sprintf('SampleFrame_Run%d', i);
        
        plot_save_mask(vorticity, b_mask, xcorners, ycorners, ...
                      setup.figures.axisFontSize, setup.figures.titleFontSize, ...
                      sample_frame_dir, setup.instantaneous.runs, ...
                      [setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)], ...
                      i, sample_variableName, 'Inst', sample_frame_dir, ...
                      vort_lower_limit, vort_upper_limit);
        
        % Load and enhance sample frame with gamma1 contours
        sample_fig_filename = fullfile(sample_frame_dir, [sample_variableName, num2str(setup.instantaneous.windowSize(i,1)), 'x', num2str(setup.instantaneous.windowSize(i,2)), '.fig']);
        fig = openfig(sample_fig_filename);
        ax = gca;
        hold(ax, 'on');
        
        [C_pos, h_pos] = contour(ax, x, y, gammaField, [0.5, 1, 2], 'LineWidth', 2, 'LineStyle', '-');
        [C_neg, h_neg] = contour(ax, x, y, gammaField, [-2, -1, -0.5], 'LineWidth', 2, 'LineStyle', '--');
        h_pos.LineColor = 'yellow';
        h_neg.LineColor = 'cyan';
        
%         title(ax, sprintf('Sample Frame: Instantaneous Vorticity with $\\Gamma_1$ Contours\\nFrame %d, Run %d', ...
%                      mid_frame, i), ...
%               'Interpreter', 'latex', 'FontSize', setup.figures.titleFontSize);
        
        xlabel(ax, 'x [mm]', 'Interpreter', 'latex', 'FontSize', setup.figures.labelFontSize);
        ylabel(ax, 'y [mm]', 'Interpreter', 'latex', 'FontSize', setup.figures.labelFontSize);
        
        cb = colorbar(ax);
        cb.Label.String = 'Vorticity [1/s]';
        cb.Label.Interpreter = 'latex';
        cb.Label.FontSize = setup.figures.labelFontSize;
        
        hold(ax, 'off');
        
        % Save enhanced sample frame with safe filename
        enhanced_base = sprintf('SampleFrame_Run%d_%dx%d_enhanced', i, setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2));
        saveas(fig, fullfile(sample_frame_dir, [enhanced_base, '.fig']));
        saveas(fig, fullfile(sample_frame_dir, [enhanced_base, '.png']));
        saveas(fig, fullfile(sample_frame_dir, [enhanced_base, '.eps']));
        close(fig);
    end
    
    fprintf('Vortex video generation completed at %s\n', datetime('now'));
end

end