function make_streamline_video(setup, endpoint, CameraNo, type, use_merged, acq_freq)
% MAKE_STREAMLINE_VIDEO - Create video of streamlines from velocity fields
%   Usage: make_streamline_video(setup, endpoint, CameraNo, type, use_merged, acq_freq)
%   acq_freq: acquisition frequency in Hz

if nargin < 4
    type = 'Instantaneous';
end
if nargin < 5
    use_merged = false;
end
if nargin < 6
    error('Acquisition frequency (acq_freq) must be provided.');
end

if ~setup.pipeline.streamline_video
    return
end

fprintf('Creating streamline video for base: %s at %s\n', setup.directory.base, datetime('now'));

paths = get_data_paths(setup, type, endpoint, CameraNo, use_merged);
dataloc = paths.data_dir;
output_dir = fullfile(paths.video_dir, 'Streamline_Analysis');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

for i = setup.instantaneous.runs
    fprintf('Processing run %d for streamline video...\n', i);

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

    % Video setup
    video_filename = fullfile(output_dir, sprintf('Streamline_Video_Run%d_%dx%d.mp4', ...
                             i, setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)));
    v = VideoWriter(video_filename, 'MPEG-4');
    v.FrameRate = 30;
    open(v);

    for imNo = 1:setup.imProperties.caseImages
        VelData = load(fullfile(dataloc, sprintf('%05d.mat', imNo)));
        u = VelData.piv_result(i).ux;
        v_vel = VelData.piv_result(i).uy;

        % Mask
        u(b_mask) = NaN;
        v_vel(b_mask) = NaN;

        % --- Preferred figure and axis setup ---
        fig = figure('Visible', 'off');
        set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
        ax1 = axes;
        % Plot streamlines
        h = streamslice(x, y, u, v_vel, 7);
        set(h, 'Color', 'k');
        title(sprintf('Streamlines Frame %d, Run %d, Time: %.3fs', imNo, i, (imNo-1)/acq_freq), ...
              'FontSize', setup.figures.titleFontSize, 'Interpreter', 'latex');
        daspect([1 1 1]);
        xlim([x(1,1) x(1,end)]);
        ylim([y(end,1) y(1,1)]);
        set(ax1, 'YDir', 'normal', 'FontSize', setup.figures.axisFontSize, ...
            'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');
        hold(ax1, 'on');

        % Mask overlay (light grey, transparent)
        mask_vis = double(b_mask);
        [rows, cols] = size(mask_vis);
        mask_rgb = cat(3, ones(rows, cols) * 0.8, ones(rows, cols) * 0.8, ones(rows, cols) * 0.8);
        j = image(ax1, xcorners, ycorners, mask_rgb);
        set(j, 'AlphaData', mask_vis * 0.7);

        daspect(ax1, [1 1 1]);
        hold(ax1, 'off');

        xlabel(ax1, 'x [mm]', 'Interpreter', 'latex', 'FontSize', setup.figures.labelFontSize);
        ylabel(ax1, 'y [mm]', 'Interpreter', 'latex', 'FontSize', setup.figures.labelFontSize);

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

fprintf('Streamline video generation completed at %s\n', datetime('now'));
