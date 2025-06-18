function multi_camera_coord_editor(setup, type, endpoint, use_merged)
    % MULTI_CAMERA_COORD_EDITOR - Comprehensive coordinate editor for multiple cameras
    % Inputs:
    %   setup - Setup structure containing camera configuration and paths
    %   type - 'Instantaneous' or 'Ensemble' data type (optional, default: 'Instantaneous')
    %   endpoint - endpoint subfolder (optional, default: '')
    %   use_merged - logical, true to use merged data, false for traditional camera data (optional, default: false)
    
    if nargin < 2
        type = 'Instantaneous';
    end
    if nargin < 3
        endpoint = '';
    end
    if nargin < 4
        use_merged = false;
    end
    
    % Validate setup structure
    if ~isfield(setup, 'imProperties') || ~isfield(setup.imProperties, 'cameraCount')
        error('setup.imProperties.cameraCount not found. Please ensure setup structure is properly configured.');
    end
    
    camera_count = setup.imProperties.cameraCount;
    
    % Determine which data type to process
    process_instantaneous =  setup.pipeline.instantaneous_cords_multi;
    process_ensemble =  setup.pipeline.ensemble_cords_multi;
    
    if ~process_instantaneous && ~process_ensemble
        return
    end
    
    % Process instantaneous data
    if process_instantaneous
        fprintf('\n=== Processing Instantaneous Data ===\n');
        for run_idx = 1:length(setup.instantaneous.runs)
            run = setup.instantaneous.runs(run_idx);
            fprintf('Processing instantaneous run %d...\n', run);
            
            % Load data for all cameras
            camera_data = load_all_camera_data(setup, camera_count, run, 'Instantaneous', endpoint, use_merged);
            
            % Create multi-camera visualization
            create_multi_camera_display(camera_data, setup, run, 'Instantaneous', camera_count);
            
            % Interactive coordinate editing
            updated_coords = interactive_coordinate_editor(camera_data, setup, run, 'Instantaneous', camera_count);
            
            % Save updated coordinates
            save_updated_coordinates(updated_coords, setup, run, 'Instantaneous', camera_count);
            
            % Show confirmation
            show_confirmation_display(updated_coords, setup, run, 'Instantaneous', camera_count);
            
            if run_idx < length(setup.instantaneous.runs)
                input('Press Enter to continue to next instantaneous run...');
            end
        end
    end
    
    % Process ensemble data
    if process_ensemble
        fprintf('\n=== Processing Ensemble Data ===\n');
        for run_idx = 1:length(setup.ensemble.runs)
            run = setup.ensemble.runs(run_idx);
            fprintf('Processing ensemble run %d...\n', run);
            
            % Load data for all cameras
            camera_data = load_all_camera_data(setup, camera_count, run, 'Ensemble', endpoint, use_merged);
            
            % Create multi-camera visualization
            create_multi_camera_display(camera_data, setup, run, 'Ensemble', camera_count);
            
            % Interactive coordinate editing
            updated_coords = interactive_coordinate_editor(camera_data, setup, run, 'Ensemble', camera_count);
            
            % Save updated coordinates
            save_updated_coordinates(updated_coords, setup, run, 'Ensemble', camera_count);
            
            % Show confirmation
            show_confirmation_display(updated_coords, setup, run, 'Ensemble', camera_count);
            
            if run_idx < length(setup.ensemble.runs)
                input('Press Enter to continue to next ensemble run...');
            end
        end
    end
    
    fprintf('\nMulti-camera coordinate editing completed successfully at %s\n', datetime('now'));
end

function camera_data = load_all_camera_data(setup, camera_count, run, data_type, endpoint, use_merged)
    % Load PIV data and coordinates for all cameras using centralized paths
    
    if nargin < 5
        endpoint = '';
    end
    if nargin < 6
        use_merged = false;
    end
    
    camera_data = struct();
    
    for camera_no = 1:camera_count
        try
            % Get data paths using centralized function
            paths = get_data_paths(setup, data_type, endpoint, camera_no, use_merged);
            dataloc = paths.data_dir;
            
            % Load velocity data (first frame)
            if strcmp(data_type, 'Instantaneous')
                vel_file = fullfile(dataloc, sprintf(setup.instantaneous.nameConvention{1}, 1));
            else % Ensemble
                vel_file = fullfile(dataloc, sprintf(setup.ensemble.nameConvention{1}, 1));
            end
            
            if ~exist(vel_file, 'file')
                warning('Velocity file not found for Camera %d: %s', camera_no, vel_file);
                continue;
            end
            
            VelData = load(vel_file);
            
            % Load coordinates
            coord_file = fullfile(dataloc, 'Co_ords.mat');
            if ~exist(coord_file, 'file')
                warning('Coordinate file not found for Camera %d: %s', camera_no, coord_file);
                continue;
            end
            
            Co_ords = load(coord_file);
            
            % Store camera data using dynamic field names to avoid structure conflicts
            camera_field = sprintf('camera_%d', camera_no);
            camera_data.(camera_field) = VelData.piv_result(run);
            camera_data.(camera_field).coordinates = Co_ords.Co_ords(run);
            camera_data.(camera_field).dataloc = dataloc;
            camera_data.(camera_field).coord_file = coord_file;
            
            % Calculate corner coordinates
            camera_data.(camera_field).corners = calculate_corner_coordinates(Co_ords.Co_ords(run));
            
            fprintf('Loaded data for Camera %d\n', camera_no);
            
        catch ME
            warning('Failed to load data for Camera %d: %s', camera_no, ME.message);
        end
    end
    
    % Check if any cameras were loaded
    if isempty(fieldnames(camera_data))
        error('No camera data could be loaded. Please check file paths and setup configuration.');
    end
end

function corners = calculate_corner_coordinates(coords)
% Calculate corner coordinates from coordinate matrices
% Note: In MATLAB matrix indexing, (1,1) is top-left, (end,end) is bottom-right
    
    corners.top_left = [coords.x(1,1), coords.y(1,1)];           % First row, first column
    corners.top_right = [coords.x(1,end), coords.y(1,end)];     % First row, last column  
    corners.bottom_left = [coords.x(end,1), coords.y(end,1)];   % Last row, first column
    corners.bottom_right = [coords.x(end,end), coords.y(end,end)]; % Last row, last column
    
    % Calculate extents for plotting (from bottom-left to top-right)
    corners.x_range = [coords.x(end,1), coords.x(1,end)];   % From bottom-left to top-right X
    corners.y_range = [coords.y(end,1), coords.y(1,end)];   % From bottom-left to top-right Y
end

function create_multi_camera_display(camera_data, ~, run, data_type, ~)
% Create subplot display showing all cameras with coordinate annotations
    
    % Determine subplot arrangement
    camera_numbers = fieldnames(camera_data);
    camera_count = length(camera_numbers);
    subplot_rows = ceil(sqrt(camera_count));
    subplot_cols = ceil(camera_count / subplot_rows);
    
    % Create figure
    fig_title = sprintf('%s Run %d - Multi-Camera Coordinate Overview', data_type, run);
    main_fig = figure('Name', fig_title, 'Position', [100, 100, 1200, 800]);
    
    camera_numbers = fieldnames(camera_data);
    for i = 1:length(camera_numbers)
        camera_no = str2double(extractAfter(camera_numbers{i}, 'camera_'));
        cam_data = camera_data.(camera_numbers{i});
        
        % Create subplot
        subplot(subplot_rows, subplot_cols, i);
        
        % Prepare data for visualization
        data_to_show = cam_data.ux;
        
        % Handle mask (different field names in different data types)
        if isfield(cam_data, 'b_mask')
            nan_mask = cam_data.b_mask;
        elseif isfield(cam_data, 'nan_mask')
            nan_mask = cam_data.nan_mask;
        else
            nan_mask = isnan(data_to_show);
        end
        
        % Create grey background for masked regions
        [rows, cols] = size(nan_mask);
        grey_bg = cat(3, ones(rows, cols) * 0.8, ones(rows, cols) * 0.8, ones(rows, cols) * 0.8);
          % Display data
        h = imagesc(cam_data.corners.x_range, cam_data.corners.y_range, data_to_show);
        set(h, 'AlphaData', ~nan_mask);
        set(gca, 'YDir', 'normal');
        hold on;
        
        % Add grey overlay for masked regions
        h_grey = image(cam_data.corners.x_range, cam_data.corners.y_range, grey_bg);
        set(h_grey, 'AlphaData', double(nan_mask) * 0.7);        % Set initial tight axis limits to reduce empty space
        actual_data_x = cam_data.coordinates.x;
        actual_data_y = cam_data.coordinates.y;
        
        % Expand limits to accommodate annotations (15% extra space)
        annotation_x_margin = (max(actual_data_x(:)) - min(actual_data_x(:))) * 0.15;
        annotation_y_margin = (max(actual_data_y(:)) - min(actual_data_y(:))) * 0.10;
        
        xlim([min(actual_data_x(:)) - annotation_x_margin, max(actual_data_x(:)) + annotation_x_margin]);
        ylim([min(actual_data_y(:)) - annotation_y_margin, max(actual_data_y(:)) + annotation_y_margin]);
        
        % Add corner annotations with expanded axis limits
        add_corner_annotations(cam_data.corners, camera_no);
        
        % Improved formatting
        title(sprintf('Camera %d', camera_no), 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.2, 0.2, 0.2]);
        xlabel('X coordinate (mm)', 'FontSize', 11, 'FontWeight', 'normal');
        ylabel('Y coordinate (mm)', 'FontSize', 11, 'FontWeight', 'normal');
        
        % Improved colorbar
        cb = colorbar;
        cb.Label.String = 'Velocity (m/s)';
        cb.Label.FontSize = 10;
        cb.FontSize = 9;
          % Grid styling
        grid on;
        set(gca, 'GridColor', [0.8, 0.8, 0.8], 'GridAlpha', 0.6, 'GridLineStyle', '-');
        set(gca, 'LineWidth', 1, 'FontSize', 9);
        
        % Set correct aspect ratio and keep tight axis limits
        daspect([1, 1, 1]);
        axis tight;
        
        hold off;
    end
      % Add overall title and instructions
    sgtitle({fig_title; ...
            'Corner coordinates shown: BL=Bottom-Left, BR=Bottom-Right, TL=Top-Left, TR=Top-Right'; ...
            'Close this window to proceed to coordinate editing'}, ...
            'FontSize', 14, 'FontWeight', 'bold');
    
    % Wait for user to close the window
    fprintf('\n=== Coordinate Overview ===\n');
    fprintf('Review the coordinate systems for all cameras.\n');
    fprintf('Close the figure window when ready to proceed to coordinate editing.\n');
    
    waitfor(main_fig);
end

function add_corner_annotations(corners, ~)
% Add text annotations for all corner coordinates with smart positioning
    
    % Get current axis limits to position text outside data area
    ax_xlim = xlim;
    ax_ylim = ylim;
    x_range = ax_xlim(2) - ax_xlim(1);
    y_range = ax_ylim(2) - ax_ylim(1);
    
    % Calculate offset distances (as percentage of axis range)
    x_offset = x_range * 0.08;  % 8% of x-range
    y_offset = y_range * 0.05;  % 5% of y-range
    
    % Enhanced annotation styling
    font_size = 11;
    bg_color = [1, 1, 1, 0.95]; % White background with high opacity
    text_props = {'FontSize', font_size, 'FontWeight', 'bold', ...
                  'BackgroundColor', bg_color, 'EdgeColor', 'black', ...
                  'HorizontalAlignment', 'center', 'LineWidth', 1.5, ...
                  'Margin', 3};
    
    % Bottom-left - position to the left and slightly down
    bl_text_x = corners.bottom_left(1) - x_offset;
    bl_text_y = corners.bottom_left(2) - y_offset;
    text(bl_text_x, bl_text_y, ...
         sprintf('BL: (%.1f, %.1f)', corners.bottom_left(1), corners.bottom_left(2)), ...
         text_props{:}, 'Color', [0.8, 0, 0], 'HorizontalAlignment', 'right');
    
    % Draw leader line
    plot([bl_text_x + x_offset*0.3, corners.bottom_left(1)], ...
         [bl_text_y + y_offset*0.3, corners.bottom_left(2)], ...
         'k-', 'LineWidth', 1, 'Color', [0.5, 0.5, 0.5]);
    plot(corners.bottom_left(1), corners.bottom_left(2), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', [0.8, 0, 0]);
    
    % Bottom-right - position to the right and slightly down
    br_text_x = corners.bottom_right(1) + x_offset;
    br_text_y = corners.bottom_right(2) - y_offset;
    text(br_text_x, br_text_y, ...
         sprintf('BR: (%.1f, %.1f)', corners.bottom_right(1), corners.bottom_right(2)), ...
         text_props{:}, 'Color', [0, 0, 0.8], 'HorizontalAlignment', 'left');
    
    % Draw leader line
    plot([br_text_x - x_offset*0.3, corners.bottom_right(1)], ...
         [br_text_y + y_offset*0.3, corners.bottom_right(2)], ...
         'k-', 'LineWidth', 1, 'Color', [0.5, 0.5, 0.5]);
    plot(corners.bottom_right(1), corners.bottom_right(2), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', [0, 0, 0.8]);
    
    % Top-left - position to the left and slightly up
    tl_text_x = corners.top_left(1) - x_offset;
    tl_text_y = corners.top_left(2) + y_offset;
    text(tl_text_x, tl_text_y, ...
         sprintf('TL: (%.1f, %.1f)', corners.top_left(1), corners.top_left(2)), ...
         text_props{:}, 'Color', [0, 0.6, 0], 'HorizontalAlignment', 'right');
    
    % Draw leader line
    plot([tl_text_x + x_offset*0.3, corners.top_left(1)], ...
         [tl_text_y - y_offset*0.3, corners.top_left(2)], ...
         'k-', 'LineWidth', 1, 'Color', [0.5, 0.5, 0.5]);
    plot(corners.top_left(1), corners.top_left(2), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', [0, 0.6, 0]);
    
    % Top-right - position to the right and slightly up
    tr_text_x = corners.top_right(1) + x_offset;
    tr_text_y = corners.top_right(2) + y_offset;
    text(tr_text_x, tr_text_y, ...
         sprintf('TR: (%.1f, %.1f)', corners.top_right(1), corners.top_right(2)), ...
         text_props{:}, 'Color', [0.7, 0, 0.7], 'HorizontalAlignment', 'left');
    
    % Draw leader line
    plot([tr_text_x - x_offset*0.3, corners.top_right(1)], ...
         [tr_text_y - y_offset*0.3, corners.top_right(2)], ...
         'k-', 'LineWidth', 1, 'Color', [0.5, 0.5, 0.5]);
    plot(corners.top_right(1), corners.top_right(2), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', [0.7, 0, 0.7]);
end

function updated_coords = interactive_coordinate_editor(camera_data, ~, ~, ~, ~)
% Interactive editor for setting bottom-left coordinates
    
    updated_coords = camera_data;
    camera_numbers = fieldnames(camera_data);
    
    fprintf('\n=== Interactive Coordinate Editor ===\n');
    fprintf('You will now set the bottom-left (origin) coordinates for each camera.\n');
    fprintf('For each camera, you can:\n');
    fprintf('  1. Enter new X,Y coordinates manually\n');
    fprintf('  2. Keep current coordinates (press Enter)\n');
    fprintf('  3. Use graphical point selection\n\n');
    
    for i = 1:length(camera_numbers)
        camera_no = str2double(camera_numbers{i});
        cam_data = camera_data.(camera_numbers{i});
        
        fprintf('--- Camera %d ---\n', camera_no);
        fprintf('Current bottom-left coordinates: (%.2f, %.2f)\n', ...
                cam_data.corners.bottom_left(1), cam_data.corners.bottom_left(2));
        
        % Get user input
        choice = input('Enter method: [1] Manual entry, [2] Keep current, [3] Graphical selection: ');
        
        new_x = cam_data.corners.bottom_left(1);
        new_y = cam_data.corners.bottom_left(2);
        
        switch choice
            case 1
                % Manual coordinate entry
                new_coords = input('Enter new bottom-left coordinates [x, y]: ');
                if ~isempty(new_coords) && length(new_coords) == 2
                    new_x = new_coords(1);
                    new_y = new_coords(2);
                end
                
            case 2
                % Keep current coordinates
                fprintf('Keeping current coordinates for Camera %d\n', camera_no);
                
            case 3
                % Graphical selection
                [new_x, new_y] = graphical_coordinate_selection(cam_data, camera_no);
                
            otherwise
                fprintf('Invalid choice. Keeping current coordinates.\n');
        end
        
        % Update coordinates if changed
        if new_x ~= cam_data.corners.bottom_left(1) || new_y ~= cam_data.corners.bottom_left(2)
            updated_coords.(camera_numbers{i}) = update_camera_coordinates(cam_data, new_x, new_y);
            fprintf('Updated Camera %d bottom-left to: (%.2f, %.2f)\n', camera_no, new_x, new_y);
        end
    end
end

function [new_x, new_y] = graphical_coordinate_selection(cam_data, camera_no)
% Graphical interface for coordinate selection
    
    % Create figure for this camera
    fig = figure('Name', sprintf('Camera %d - Click to set bottom-left origin', camera_no), ...
                 'Position', [200, 200, 800, 600]);
    
    % Prepare visualization data
    data_to_show = cam_data.ux;
    
    % Handle mask
    if isfield(cam_data, 'b_mask')
        nan_mask = cam_data.b_mask;
    elseif isfield(cam_data, 'nan_mask')
        nan_mask = cam_data.nan_mask;
    else
        nan_mask = isnan(data_to_show);
    end
    
    % Display data
    h = imagesc(cam_data.corners.x_range, cam_data.corners.y_range, data_to_show);
    set(h, 'AlphaData', ~nan_mask);
    set(gca, 'YDir', 'normal');
    hold on;
      % Add grey overlay for masked regions
    [rows, cols] = size(nan_mask);
    grey_bg = cat(3, ones(rows, cols) * 0.8, ones(rows, cols) * 0.8, ones(rows, cols) * 0.8);
    h_grey = image(cam_data.corners.x_range, cam_data.corners.y_range, grey_bg);
    set(h_grey, 'AlphaData', double(nan_mask) * 0.7);
      % Set expanded axis limits to accommodate annotations
    actual_data_x = cam_data.coordinates.x;
    actual_data_y = cam_data.coordinates.y;
    annotation_x_margin = (max(actual_data_x(:)) - min(actual_data_x(:))) * 0.15;
    annotation_y_margin = (max(actual_data_y(:)) - min(actual_data_y(:))) * 0.10;
    xlim([min(actual_data_x(:)) - annotation_x_margin, max(actual_data_x(:)) + annotation_x_margin]);
    ylim([min(actual_data_y(:)) - annotation_y_margin, max(actual_data_y(:)) + annotation_y_margin]);
    
    % Add current corner annotations
    add_corner_annotations(cam_data.corners, camera_no);
    
    % Add crosshairs at current origin
    plot(cam_data.corners.bottom_left(1), cam_data.corners.bottom_left(2), ...
         'r+', 'MarkerSize', 15, 'LineWidth', 3);
    
    % Improved formatting
    title(sprintf('Camera %d - Use zoom/pan tools, then click to set new bottom-left origin', camera_no), ...
          'FontSize', 13, 'FontWeight', 'bold', 'Color', [0.2, 0.2, 0.2]);
    xlabel('X coordinate (mm)', 'FontSize', 11, 'FontWeight', 'normal');
    ylabel('Y coordinate (mm)', 'FontSize', 11, 'FontWeight', 'normal');
    
    % Improved colorbar
    cb = colorbar;
    cb.Label.String = 'Velocity (m/s)';
    cb.Label.FontSize = 10;
      % Grid styling
    grid on;
    set(gca, 'GridColor', [0.8, 0.8, 0.8], 'GridAlpha', 0.6);
    
    % Set correct aspect ratio
    daspect([1, 1, 1]);
    
    % Instructions
    fprintf('\nGraphical selection for Camera %d:\n', camera_no);
    fprintf('1. Use zoom and pan tools to navigate\n');
    fprintf('2. Press ENTER in command window when ready\n');
    fprintf('3. Click once on the plot to set new origin\n');
    
    % Wait for user to be ready
    input('Press ENTER when ready to select point, then click on plot: ');
    
    % Get click coordinates
    [new_x, new_y] = ginput(1);
    
    % Close figure
    close(fig);
    
    if isempty(new_x) || isempty(new_y)
        % User cancelled, keep original coordinates
        new_x = cam_data.corners.bottom_left(1);
        new_y = cam_data.corners.bottom_left(2);
        fprintf('Selection cancelled. Keeping original coordinates.\n');
    else
        fprintf('Selected new bottom-left coordinates: (%.2f, %.2f)\n', new_x, new_y);
    end
end

function updated_cam_data = update_camera_coordinates(cam_data, new_bottom_left_x, new_bottom_left_y)
% Update coordinate system based on new bottom-left coordinates
    
    % Calculate offset
    current_x = cam_data.corners.bottom_left(1);
    current_y = cam_data.corners.bottom_left(2);
    
    x_offset = current_x - new_bottom_left_x;
    y_offset = current_y - new_bottom_left_y;
    
    % Apply offset to coordinate matrices
    updated_cam_data = cam_data;
    updated_cam_data.coordinates.x = cam_data.coordinates.x - x_offset;
    updated_cam_data.coordinates.y = cam_data.coordinates.y - y_offset;
    
    % Update corner coordinates
    updated_cam_data.corners = calculate_corner_coordinates(updated_cam_data.coordinates);
end

function save_updated_coordinates(updated_coords, ~, run, ~, ~)
% Save updated coordinates to files
    
    camera_numbers = fieldnames(updated_coords);
    
    fprintf('\nSaving updated coordinates...\n');
    for i = 1:length(camera_numbers)
        camera_no = str2double(extractAfter(camera_numbers{i}, 'camera_'));
        cam_data = updated_coords.(camera_numbers{i});
        
        try
            % Load existing coordinate file
            Co_ords = load(cam_data.coord_file);
            
            % Update the specific run
            Co_ords.Co_ords(run) = cam_data.coordinates;
            
            % Save back to file
            save(cam_data.coord_file, '-struct', 'Co_ords');
            
            fprintf('Saved coordinates for Camera %d\n', camera_no);
            
        catch ME
            warning('Failed to save coordinates for Camera %d: %s', camera_no, ME.message);
        end
    end
end

function show_confirmation_display(updated_coords, ~, run, data_type, ~)
% Show confirmation display with updated coordinates
    
    % Determine subplot arrangement
    camera_numbers = fieldnames(updated_coords);
    camera_count = length(camera_numbers);
    subplot_rows = ceil(sqrt(camera_count));
    subplot_cols = ceil(camera_count / subplot_rows);
    
    % Create confirmation figure
    fig_title = sprintf('%s Run %d - Updated Coordinates Confirmation', data_type, run);
    conf_fig = figure('Name', fig_title, 'Position', [150, 150, 1200, 800]);
    
    camera_numbers = fieldnames(updated_coords);
    
    for i = 1:length(camera_numbers)
        camera_no = str2double(camera_numbers{i});
        cam_data = updated_coords.(camera_numbers{i});
        
        % Create subplot
        subplot(subplot_rows, subplot_cols, i);
        
        % Prepare data for visualization
        data_to_show = cam_data.ux;
        
        % Handle mask
        if isfield(cam_data, 'b_mask')
            nan_mask = cam_data.b_mask;
        elseif isfield(cam_data, 'nan_mask')
            nan_mask = cam_data.nan_mask;
        else
            nan_mask = isnan(data_to_show);
        end
        
        % Display data with updated coordinates
        h = imagesc(cam_data.corners.x_range, cam_data.corners.y_range, data_to_show);
        set(h, 'AlphaData', ~nan_mask);
        set(gca, 'YDir', 'normal');
        hold on;
          % Add grey overlay for masked regions
        [rows, cols] = size(nan_mask);
        grey_bg = cat(3, ones(rows, cols) * 0.8, ones(rows, cols) * 0.8, ones(rows, cols) * 0.8);
        h_grey = image(cam_data.corners.x_range, cam_data.corners.y_range, grey_bg);
        set(h_grey, 'AlphaData', double(nan_mask) * 0.7);
          % Set expanded axis limits to accommodate annotations
        actual_data_x = cam_data.coordinates.x;
        actual_data_y = cam_data.coordinates.y;
        annotation_x_margin = (max(actual_data_x(:)) - min(actual_data_x(:))) * 0.15;
        annotation_y_margin = (max(actual_data_y(:)) - min(actual_data_y(:))) * 0.10;
        xlim([min(actual_data_x(:)) - annotation_x_margin, max(actual_data_x(:)) + annotation_x_margin]);
        ylim([min(actual_data_y(:)) - annotation_y_margin, max(actual_data_y(:)) + annotation_y_margin]);
        
        % Add updated corner annotations
        add_corner_annotations(cam_data.corners, camera_no);
        
        % Add origin crosshairs
        plot(0, 0, 'r+', 'MarkerSize', 15, 'LineWidth', 3);
        plot([min(actual_data_x(:)), max(actual_data_x(:))], [0, 0], 'r--', 'LineWidth', 1, 'Color', [0.8, 0, 0]);
        plot([0, 0], [min(actual_data_y(:)), max(actual_data_y(:))], 'r--', 'LineWidth', 1, 'Color', [0.8, 0, 0]);
        
        % Improved formatting
        title(sprintf('Camera %d - UPDATED', camera_no), 'FontSize', 14, 'FontWeight', 'bold', 'Color', [0, 0.6, 0]);
        xlabel('X coordinate (mm)', 'FontSize', 11, 'FontWeight', 'normal');
        ylabel('Y coordinate (mm)', 'FontSize', 11, 'FontWeight', 'normal');
        
        % Improved colorbar
        cb = colorbar;
        cb.Label.String = 'Velocity (m/s)';
        cb.Label.FontSize = 10;
        cb.FontSize = 9;
          % Grid styling
        grid on;
        set(gca, 'GridColor', [0.8, 0.8, 0.8], 'GridAlpha', 0.6);
        
        % Set correct aspect ratio and keep tight axis limits
        daspect([1, 1, 1]);
        axis tight;
        
        hold off;
    end
    
    % Add overall title
    sgtitle({fig_title; ...
            'Red crosshairs show new coordinate origins'; ...
            'This window will close automatically in 5 seconds'}, ...
            'FontSize', 14, 'FontWeight', 'bold');
    
    % Display for 5 seconds then close
    pause(5);
    close(conf_fig);
    
    fprintf('Coordinate update confirmation displayed.\n');
end
