function merge_hamming(setup, type)
    % MERGE_HAMMING - Merges calibrated vector fields from multiple cameras
    % Inputs:
    %   setup - setup structure containing camera and processing parameters
    %   type - 'Calibrated' for instantaneous data, 'SumCalibrated' for ensemble data
    
    % Validate inputs and extract parameters
    params = validate_and_extract_params(setup, type);
    
    % Load coordinate systems and calculate overlap regions
    [Co_ords, overlap_info] = setup_coordinate_systems(params);
    
    % Process each image
    process_images(params, Co_ords, overlap_info);
    
    % Save coordinate system
    save(fullfile(params.output_dir, 'Co_ords.mat'), "Co_ords");
    
    fprintf('Merge completed successfully at %s\n', datetime('now'));
end

function params = validate_and_extract_params(setup, type)
    % Extract and validate parameters from setup structure
    
    % Input validation
    if ~ischar(type) && ~isstring(type)
        error('Type must be a string: ''Calibrated'' or ''SumCalibrated''');
    end
    
    % Extract parameters
    params.type = type;
    params.camera_count = setup.camera.count;
    params.runs_merge = setup.instantaneous.runs;
    params.N = setup.imProperties.imageCount;
    params.calibrated_piv_dir = fullfile(setup.directory.base, 'CalibratedPIV');
    params.output_dir = fullfile(setup.directory.base, 'Merged');
    params.wsize = setup.camera.wsize;
    
    % Determine image count and paths based on type
    if strcmp(type, 'SumCalibrated')
        params.N_count = 1;
        params.is_ensemble = true;
        params.data_subdir = 'Ensemble';
    elseif strcmp(type, 'Calibrated')
        params.N_count = params.N;
        params.is_ensemble = false;
        params.data_subdir = 'Instantaneous';
    else
        error('Unknown type: %s. Must be ''Calibrated'' or ''SumCalibrated''', type);
    end
    
    % Create output directory
    if ~exist(params.output_dir, 'dir')
        mkdir(params.output_dir);
    end
    
    % Set truncation parameters based on window size
    params.truncation = calculate_truncation_params(params.wsize, params.runs_merge);
end

function truncation = calculate_truncation_params(wsize, runs_merge)
    % Calculate truncation parameters based on window size
    truncation = struct();
    
    for i = runs_merge
        if any(wsize(i,:) < 64)
            truncation(i).top = 5;
            truncation(i).bottom = 5;
            truncation(i).side = 5;
        else
            truncation(i).top = 0;
            truncation(i).bottom = 0;
            truncation(i).side = 0;
        end
    end
end

function [Co_ords, overlap_info] = setup_coordinate_systems(params)
    % Load coordinate systems and calculate overlap regions
    
    Co_ords = struct();
    overlap_info = struct();
    pos = struct();
    
    % Load coordinate data for each run and camera
    for i = params.runs_merge
        pos(i).x = cell(1, params.camera_count);
        pos(i).y = cell(1, params.camera_count);
        
        for camera_no = 1:params.camera_count
            coords_path = fullfile(params.calibrated_piv_dir, num2str(params.N), ...
                                 ['Cam', num2str(camera_no)], params.data_subdir, 'Co_ords.mat');
            
            if ~exist(coords_path, 'file')
                error('Coordinate file not found: %s', coords_path);
            end
            
            data = load(coords_path);
            
            % Apply truncation
            trunc = params.truncation(i);
            pos(i).x{camera_no} = data.Co_ords(i).x(trunc.top+1:end-trunc.bottom, trunc.side+1:end-trunc.side);
            pos(i).y{camera_no} = data.Co_ords(i).y(trunc.top+1:end-trunc.bottom, trunc.side+1:end-trunc.side);
        end
        
        % Calculate overlap regions and Hanning windows
        overlap_info(i) = calculate_overlap_regions(pos(i), params.camera_count);
        
        % Create merged coordinate system
        Co_ords(i) = create_merged_coordinates(pos(i), params.camera_count);
        
        % Store camera coordinates for interpolation
        Co_ords(i).camera_coords = pos(i);
    end
end

function overlap_info = calculate_overlap_regions(pos, camera_count)
    % Calculate overlap regions between adjacent cameras
    
    overlap_info.regions = zeros(1, camera_count-1);
    overlap_info.hanning_windows = cell(1, camera_count-1);
    
    for camera_no = 1:camera_count-1
        if camera_no < camera_count
            % Find overlap between current and next camera
            [~, overlap_idx] = min(abs(pos.x{camera_no}(1,:) - pos.x{camera_no+1}(1,1)));
            overlap_info.regions(camera_no) = length(pos.x{camera_no}) - overlap_idx;
            
            % Create Hanning window for smooth transition
            if overlap_info.regions(camera_no) > 0
                h = hanning(overlap_info.regions(camera_no) * 2 + 1);
                overlap_info.hanning_windows{camera_no} = h(1:overlap_info.regions(camera_no))';
            else
                overlap_info.hanning_windows{camera_no} = [];
                warning('No overlap found between cameras %d and %d', camera_no, camera_no+1);
            end
        end
    end
end

function Co_ords = create_merged_coordinates(pos, camera_count)
    % Create merged coordinate system spanning all cameras
    
    % Concatenate all coordinates
    x_merge = cat(2, pos.x{:});
    y_merge = cat(2, pos.y{:});
    
    % Find domain bounds
    bounds = calculate_domain_bounds(x_merge, y_merge);
    
    % Calculate cell density
    cell_density_x = min(min(diff(pos.x{1}, 1, 2)));
    cell_density_y = -min(min(diff(pos.y{1}, 1, 1)));
    
    % Create merged grid
    x_lims = [bounds.min_x, bounds.max_x];
    y_lims = [bounds.min_y, bounds.max_y];
    
    [x, y] = meshgrid(x_lims(1):cell_density_x:x_lims(2), ...
                      y_lims(1):cell_density_y:y_lims(2));
    
    Co_ords.x = double(x);
    Co_ords.y = flipud(double(y));
end

function bounds = calculate_domain_bounds(x_merge, y_merge)
    % Calculate domain bounds from merged coordinates
    
    % Find bounds that minimize distance from zero (closest to centerline)
    top_row = x_merge(1, :);
    [~, pos_top] = min(abs(top_row));
    bounds.max_y = top_row(pos_top);
    
    right_col = x_merge(:, end);
    [~, pos_right] = min(abs(right_col));
    bounds.max_x = right_col(pos_right);
    
    left_col = x_merge(:, 1);
    [~, pos_left] = min(abs(left_col));
    bounds.min_x = left_col(pos_left);
    
    bottom_row = y_merge(end, :);
    [~, pos_bottom] = max(abs(bottom_row));
    bounds.min_y = bottom_row(pos_bottom);
end

function process_images(params, Co_ords, overlap_info)
    % Process all images for merging
    
    parfor im_no = 1:params.N_count
        fprintf('Processing image %d of %d at %s\n', im_no, params.N_count, datetime('now'));
        
        merged_results = struct();
        
        for i = params.runs_merge
            merged_results(i) = process_single_image(im_no, i, params, Co_ords(i), overlap_info(i));
        end
        
        % Save merged results
        output_file = fullfile(params.output_dir, sprintf('%05d.mat', im_no));
        save_piv_results(output_file, merged_results);
    end
end

function piv_result = process_single_image(im_no, run_idx, params, coords, overlap_info)
    % Process a single image for a specific run
    
    % Initialize result structure
    piv_result = initialize_piv_result(coords, params.is_ensemble);
    
    % Load and process data from each camera
    camera_data = load_camera_data(im_no, run_idx, params);
    
    % Apply Hanning window weighting in overlap regions
    weighted_data = apply_hanning_weighting(camera_data, overlap_info, params.camera_count);
    
    % Interpolate data onto merged grid
    interpolated_data = interpolate_to_merged_grid(weighted_data, coords, params.camera_count);
    
    % Sum contributions from all cameras
    piv_result = sum_camera_contributions(piv_result, interpolated_data, params.is_ensemble);
    
    % Post-process masks
    piv_result = post_process_masks(piv_result, params.is_ensemble);
end

function piv_result = initialize_piv_result(coords, is_ensemble)
    % Initialize PIV result structure
    
    grid_size = size(coords.x);
    piv_result.ux = zeros(grid_size);
    piv_result.uy = zeros(grid_size);
    piv_result.b_mask = zeros(grid_size);
    
    if is_ensemble
        piv_result.Uturb = zeros(grid_size);
        piv_result.Vturb = zeros(grid_size);
        piv_result.UturbVturb = zeros(grid_size);
        piv_result.NanMask = zeros(grid_size);
    end
end

function camera_data = load_camera_data(im_no, run_idx, params)
    % Load velocity data from all cameras
    
    camera_data = struct();
    
    for camera_no = 1:params.camera_count
        vel_path = fullfile(params.calibrated_piv_dir, num2str(params.N), ...
                           ['Cam', num2str(camera_no)], params.data_subdir, ...
                           sprintf('%05d.mat', im_no));
        
        if ~exist(vel_path, 'file')
            error('Velocity file not found: %s', vel_path);
        end
        
        vel_data = load(vel_path);
        trunc = params.truncation(run_idx);
        
        % Extract and truncate data
        camera_data(camera_no) = extract_camera_data(vel_data.piv_result(run_idx), trunc, params.is_ensemble);
    end
end

function data = extract_camera_data(piv_result, trunc, is_ensemble)
    % Extract and truncate data from PIV result
    
    data.u = piv_result.ux(trunc.top+1:end-trunc.bottom, trunc.side+1:end-trunc.side);
    data.v = piv_result.uy(trunc.top+1:end-trunc.bottom, trunc.side+1:end-trunc.side);
    data.b_mask = double(piv_result.b_mask(trunc.top+1:end-trunc.bottom, trunc.side+1:end-trunc.side));
    
    if is_ensemble
        data.Uturb = piv_result.Uturb(trunc.top+1:end-trunc.bottom, trunc.side+1:end-trunc.side);
        data.Vturb = piv_result.Vturb(trunc.top+1:end-trunc.bottom, trunc.side+1:end-trunc.side);
        data.UturbVturb = piv_result.UturbVturb(trunc.top+1:end-trunc.bottom, trunc.side+1:end-trunc.side);
        data.NanMask = double(piv_result.NanMask(trunc.top+1:end-trunc.bottom, trunc.side+1:end-trunc.side));
    end
end

function weighted_data = apply_hanning_weighting(camera_data, overlap_info, camera_count)
    % Apply Hanning window weighting in overlap regions
    
    weighted_data = camera_data;  % Copy original data
    
    for camera_no = 1:camera_count
        if camera_no > 1
            % Apply weighting to left overlap region
            overlap_size = overlap_info.regions(camera_no-1);
            if overlap_size > 0 && ~isempty(overlap_info.hanning_windows{camera_no-1})
                window = overlap_info.hanning_windows{camera_no-1};
                weighted_data(camera_no) = apply_window_to_data(weighted_data(camera_no), window, 'left', overlap_size);
            end
        end
        
        if camera_no < camera_count
            % Apply weighting to right overlap region
            overlap_size = overlap_info.regions(camera_no);
            if overlap_size > 0 && ~isempty(overlap_info.hanning_windows{camera_no})
                window = 1 - overlap_info.hanning_windows{camera_no};
                weighted_data(camera_no) = apply_window_to_data(weighted_data(camera_no), window, 'right', overlap_size);
            end
        end
    end
end

function data = apply_window_to_data(data, window, side, overlap_size)
    % Apply windowing to data fields
    
    fields = fieldnames(data);
    for f = 1:length(fields)
        field_name = fields{f};
        if strcmp(side, 'left')
            data.(field_name)(:, 1:overlap_size) = data.(field_name)(:, 1:overlap_size) .* window;
        else % right
            data.(field_name)(:, end-overlap_size+1:end) = data.(field_name)(:, end-overlap_size+1:end) .* window;
        end
    end
end

function interpolated_data = interpolate_to_merged_grid(weighted_data, coords, camera_count)
    % Interpolate camera data to merged grid
    
    interpolated_data = struct();
    
    for camera_no = 1:camera_count
        data = weighted_data(camera_no);
        camera_coords = coords.camera_coords;
        
        % Interpolate each field using scatteredInterpolant
        interpolated_data(camera_no) = interpolate_camera_data(data, camera_coords.x{camera_no}, ...
                                                              camera_coords.y{camera_no}, coords);
    end
end

function interp_data = interpolate_camera_data(data, cam_x, cam_y, merged_coords)
    % Interpolate single camera data to merged grid
    
    fields = fieldnames(data);
    interp_data = struct();
    
    for f = 1:length(fields)
        field_name = fields{f};
        field_data = data.(field_name);
        
        % Create interpolant
        interp = scatteredInterpolant(cam_x(:), cam_y(:), double(field_data(:)), 'natural', 'none');
        
        % Interpolate to merged grid
        interp_data.(field_name) = interp(merged_coords.x, merged_coords.y);
        
        % Replace NaNs with zeros
        nan_indices = isnan(interp_data.(field_name));
        interp_data.(field_name)(nan_indices) = 0;
    end
end

function piv_result = sum_camera_contributions(piv_result, interpolated_data, is_ensemble)
    % Sum contributions from all cameras
    
    camera_numbers = fieldnames(interpolated_data);
    
    for c = 1:length(camera_numbers)
        data = interpolated_data.(camera_numbers{c});
        
        piv_result.ux = piv_result.ux + data.u;
        piv_result.uy = piv_result.uy + data.v;
        piv_result.b_mask = piv_result.b_mask + data.b_mask;
        
        if is_ensemble
            piv_result.Uturb = piv_result.Uturb + data.Uturb;
            piv_result.Vturb = piv_result.Vturb + data.Vturb;
            piv_result.UturbVturb = piv_result.UturbVturb + data.UturbVturb;
            piv_result.NanMask = piv_result.NanMask + data.NanMask;
        end
    end
end

function piv_result = post_process_masks(piv_result, is_ensemble)
    % Post-process binary masks
    
    piv_result.b_mask = piv_result.b_mask > 0.01;
    
    if is_ensemble
        piv_result.NanMask = piv_result.NanMask > 0.01;
    end
end

function save_piv_results(filename, piv_results)
    % Save PIV results to file
    
    piv_result = piv_results;  % Rename for compatibility
    save(filename, 'piv_result', '-v7.3');
end



