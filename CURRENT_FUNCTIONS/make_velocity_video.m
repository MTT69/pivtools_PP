function make_velocity_video(setup, endpoint, CameraNo, type, use_merged, component, acq_freq)
% MAKE_VELOCITY_VIDEO - Create video for ux or uy velocity fields
%   Usage: make_velocity_video(setup, endpoint, CameraNo, type, use_merged, component, acq_freq)
%   component: 'ux' or 'uy'
%   acq_freq: acquisition frequency in Hz

if nargin < 4
    type = 'Instantaneous';
end
if nargin < 5
    use_merged = false;
end
if nargin < 6
    component = 'ux'; % default to ux
end
if nargin < 7
    error('Acquisition frequency (acq_freq) must be provided.');
end
if ~setup.pipeline.velocity_video
    return
end
fprintf('Creating %s velocity video for base: %s at %s\n', component, setup.directory.base, datetime('now'));

paths = get_data_paths(setup, type, endpoint, CameraNo, use_merged);
dataloc = paths.data_dir;
output_dir = fullfile(paths.video_dir, 'Velocity_Analysis', component);

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

for i = setup.instantaneous.runs
    fprintf('Processing run %d for %s velocity video...\n', i, component);

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

    % Pre-calculate limits for consistent colorbar scaling
    fprintf('Calculating %s limits for consistent scaling...\n', component);
    sample_frames = 1:50;
    all_values = [];

    for sample_idx = 1:length(sample_frames)
        sample_imNo = sample_frames(sample_idx);
        VelData_sample = load(fullfile(dataloc, sprintf('%05d.mat', sample_imNo)));
        u_sample = VelData_sample.piv_result(i).ux;
        v_sample = VelData_sample.piv_result(i).uy;

        if strcmp(component, 'ux')
            field_sample = u_sample;
        else
            field_sample = v_sample;
        end


        valid_vals = field_sample(~isnan(field_sample));
        all_values = [all_values; valid_vals(:)];

        if mod(sample_idx, 10) == 0
            fprintf('Sampled %d/%d frames for limit calculation\n', sample_idx, length(sample_frames));
        end
    end

    % Calculate consistent limits using percentiles
    lower_limit = prctile(all_values, 5);
    upper_limit = prctile(all_values, 95);

    % Ensure symmetric limits around zero
    max_abs_limit = max(abs(lower_limit), abs(upper_limit));
    lower_limit = -max_abs_limit;
    upper_limit = max_abs_limit;

    lower_limit = -0.025;
    upper_limit = 0.025;

    fprintf('%s limits: [%.2f, %.2f]\n', component, lower_limit, upper_limit);

    % Video setup
    video_filename = fullfile(output_dir, sprintf('%s_Video_Run%d_%dx%d.mp4', ...
                             component, i, setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)));
    v = VideoWriter(video_filename, 'MPEG-4');
    v.FrameRate = 30;
    open(v);

    for imNo = 1:setup.imProperties.caseImages
        VelData = load(fullfile(dataloc, sprintf('%05d.mat', imNo)));
        u = VelData.piv_result(i).ux;
        v_vel = VelData.piv_result(i).uy;

        if strcmp(component, 'ux')
            field = u;
        else
            field = v_vel;
        end

        field(b_mask) = 0;

        % Create figure using plot_save_mask with return_figure option
        variableName = sprintf('%s Velocity Field', component);
        fig = plot_save_mask(field, b_mask, xcorners, ycorners, ...
                  setup.figures.axisFontSize, setup.figures.titleFontSize, ...
                  output_dir, setup.instantaneous.runs, ...
                  [setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)], ...
                  i, variableName, 'Inst', output_dir, ...
                  lower_limit, upper_limit, true);

        set(fig, 'Visible', 'off');
        set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
        ax = gca;
        % Add time to title
        time_s = (imNo - 1) / acq_freq;
        title(ax, sprintf('%s Velocity Field, Time: %.3fs', component, time_s), ...
              'Interpreter', 'latex', 'FontSize', setup.figures.titleFontSize);

        xlabel(ax, 'x [mm]', 'Interpreter', 'latex', 'FontSize', setup.figures.labelFontSize);
        ylabel(ax, 'y [mm]', 'Interpreter', 'latex', 'FontSize', setup.figures.labelFontSize);

        cb = colorbar(ax);
        cb.Label.String = sprintf('%s velocity', component);
        cb.Label.Interpreter = 'latex';
        cb.Label.FontSize = setup.figures.labelFontSize;

        frame = getframe(fig);
        writeVideo(v, frame);
        close(fig);

        if mod(imNo, 10) == 0
            fprintf('Processed frame %d/%d\n', imNo, setup.imProperties.caseImages);
        end
    end

    close(v);
    fprintf('Video saved: %s\n', video_filename);

   
end

fprintf('%s velocity video generation completed at %s\n', component, datetime('now'));
end
