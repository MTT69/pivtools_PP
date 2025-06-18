function merge(setup, type, endpoint)
    % MERGE_HAMMING - Merges calibrated vector fields from multiple cameras
    % Inputs:
    %   setup - setup structure containing camera and processing parameters
    %   type - 'Instantaneous' for instantaneous data, 'Ensemble' for ensemble data
    %   endpoint - endpoint subfolder (optional, defaults to '')
    
    if nargin < 3
        endpoint = ''; % Default endpoint
    end
    
    % Validate inputs and extract parameters
    if setup.pipeline.merge
        
        params = validate_and_extract_params(setup, type, endpoint);
        
        % Load coordinate systems and calculate overlap regions
        [Co_ords, overlap_info] = setup_coordinate_systems(params);
        
        % Process each image
        process_images(params, Co_ords, overlap_info);
        
        % Save coordinate system
        save(fullfile(params.output_dir, 'Co_ords.mat'), "Co_ords");
        
        fprintf('Merge completed successfully at %s\n', datetime('now'));
    end
end

function params = validate_and_extract_params(setup, type, endpoint)
    % Extract and validate parameters from setup structure
    
    % Input validation
    if ~ischar(type) && ~isstring(type)
        error('Type must be a string: ''Instantaneous'' or ''Ensemble''');
    end
    
    % Extract parameters
    params.type = type;
    params.endpoint = endpoint;
    params.camera_count = setup.imProperties.cameraCount;
    params.runs_merge = setup.instantaneous.runs;
    params.N = setup.imProperties.imageCount;
    params.calibrated_piv_dir = fullfile(setup.directory.base, 'CalibratedPIV');
    
    % Create output directory path with type and endpoint
    if isempty(endpoint)
        params.output_dir = fullfile(setup.directory.base, 'Merged', type);
    else
        params.output_dir = fullfile(setup.directory.base, 'Merged', type, endpoint);
    end
    
    params.wsize = setup.instantaneous.windowSize;
    
    % Determine image count and paths based on type
    if strcmp(type, 'Ensemble')
        params.N_count = 1;
        params.is_ensemble = true;
        params.data_subdir = 'Ensemble';
    elseif strcmp(type, 'Instantaneous')
        params.N_count = params.N;
        params.is_ensemble = false;
        params.data_subdir = 'Instantaneous';
    else
        error('Unknown type: %s. Must be ''Instantaneous'' or ''Ensemble''', type);
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
            truncation(i).top = 2;
            truncation(i).bottom = 2;
            truncation(i).side = 2;
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
            % Build coordinate path with endpoint support
            if isempty(params.endpoint)
                coords_path = fullfile(params.calibrated_piv_dir, num2str(params.N), ...
                                     ['Cam', num2str(camera_no)], params.data_subdir, 'Co_ords.mat');
            else
                coords_path = fullfile(params.calibrated_piv_dir, num2str(params.N), ...
                                     ['Cam', num2str(camera_no)], params.data_subdir, params.endpoint, 'Co_ords.mat');
            end
            
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
        overlap_result = calculate_overlap_regions(pos(i), params.camera_count);
        overlap_info(i).regions = overlap_result.regions;
        overlap_info(i).hanning_windows = overlap_result.hanning_windows;
        
        % Create merged coordinate system
        coords_result = create_merged_coordinates(pos(i));
        Co_ords(i).x = coords_result.x;
        Co_ords(i).y = coords_result.y;
        Co_ords(i).camera_coords = pos(i);
    end
end

function overlap_info = calculate_overlap_regions(pos, camera_count)
    % Calculate overlap regions between adjacent cameras
    
    overlap_info.regions = zeros(1, camera_count);
    overlap_info.hanning_windows = cell(1, camera_count);
    
    for camera_no = 1:camera_count-1
        % Find overlap between current and next camera
        [~, overlap_idx] = min(abs(pos.x{camera_no}(1,:) - pos.x{camera_no+1}(1,1)));
        overlap_size = length(pos.x{camera_no}(1,:)) - overlap_idx + 1;
        overlap_info.regions(camera_no) = overlap_size;
        
        % Create Hanning window for smooth transition
        if overlap_size > 0
            h = hanning(overlap_size * 2 + 1);
            overlap_info.hanning_windows{camera_no} = h(1:overlap_size)';
        else
            overlap_info.hanning_windows{camera_no} = [];
            warning('No overlap found between cameras %d and %d', camera_no, camera_no+1);
        end
    end
end

function Co_ords = create_merged_coordinates(pos)
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
    top_row_y = y_merge(1,:);
    bounds.max_y = min(top_row_y);
    
    right_col = x_merge(:, end);
    bounds.max_x = min(right_col);
    
    left_col = x_merge(:, 1);
    bounds.min_x = max(left_col);
    
    bottom_row = y_merge(end, :);
    bounds.min_y = max(bottom_row);
end

function process_images(params, Co_ords, overlap_info)
    % Process all images for merging
    
    parfor im_no = 1:params.N_count
        
        % Create a template structure with all required fields
        template_struct = struct('ux', [], 'uy', [], 'b_mask', []);
        if params.is_ensemble
            template_struct.Uturb = [];
            template_struct.Vturb = [];
            template_struct.UturbVturb = [];
            template_struct.NanMask = [];
        end
        
        % Pre-allocate piv_result structure array with template
        max_run_idx = max(params.runs_merge);
        piv_result = repmat(template_struct, 1, max_run_idx);
        
        % Process each run and collect results
        for idx = 1:length(params.runs_merge)
            i = params.runs_merge(idx);
            piv_result(i) = process_single_image(im_no, i, params, Co_ords(i), overlap_info(i));
        end
        
        % Save merged results
        output_file = fullfile(params.output_dir, sprintf('%05d.mat', im_no));
        PIVSAVE(output_file, piv_result);
    end
end

function piv_result = process_single_image(im_no, run_idx, params, coords, overlap_info)
    % Process a single image for a specific run - collect all data then interpolate once
    
    % Initialize result structure
    piv_result = initialize_piv_result(coords, params.is_ensemble);
    
    % Load data from all cameras
    camera_data = load_camera_data(im_no, run_idx, params);
    
    % Collect all camera data into single arrays
    [all_coords, all_data] = collect_all_camera_data(camera_data, coords.camera_coords, overlap_info, params);
    
    % Single interpolation for all data
    piv_result = interpolate_all_data(all_coords, all_data, coords, params.is_ensemble);
    
    % Post-process masks
    piv_result = post_process_masks(piv_result, params.is_ensemble);
end

function [all_coords, all_data] = collect_all_camera_data(camera_data, camera_coords, overlap_info, params)
    % Collect all camera data with proper Hanning weighting
    
    all_x = [];
    all_y = [];
    all_u = [];
    all_v = [];
    all_b_mask = [];
    
    if params.is_ensemble
        all_Uturb = [];
        all_Vturb = [];
        all_UturbVturb = [];
        all_NanMask = [];
    end
    
    for camera_no = 1:params.camera_count
        % Get camera coordinates and data
        cam_x = camera_coords.x{camera_no};
        cam_y = camera_coords.y{camera_no};
        cam_data = camera_data(camera_no);
        
        % Apply Hanning weighting based on position in overlap regions
        % this method needs fixing
        % the idea is that values more central are more trusted but requried a double interpolatino i think.... just left for now uniform weights used
        % weights = calculate_hanning_weights(cam_x, camera_no, overlap_info, params.camera_count);
        weights = ones(size(cam_x));
        % Apply weights to data
        weighted_u = cam_data.u .* weights;
        weighted_v = cam_data.v .* weights;
        weighted_b_mask = cam_data.b_mask .* weights;
        
        % Collect coordinates and weighted data
        all_x = [all_x; cam_x(:)];
        all_y = [all_y; cam_y(:)];
        all_u = [all_u; weighted_u(:)];
        all_v = [all_v; weighted_v(:)];
        all_b_mask = [all_b_mask; weighted_b_mask(:)];
        
        if params.is_ensemble
            weighted_Uturb = cam_data.Uturb .* weights;
            weighted_Vturb = cam_data.Vturb .* weights;
            weighted_UturbVturb = cam_data.UturbVturb .* weights;
            weighted_NanMask = cam_data.NanMask .* weights;
            
            all_Uturb = [all_Uturb; weighted_Uturb(:)];
            all_Vturb = [all_Vturb; weighted_Vturb(:)];
            all_UturbVturb = [all_UturbVturb; weighted_UturbVturb(:)];
            all_NanMask = [all_NanMask; weighted_NanMask(:)];
        end
    end
    
    % Package results
    all_coords.x = all_x;
    all_coords.y = all_y;
    
    all_data.u = all_u;
    all_data.v = all_v;
    all_data.b_mask = all_b_mask;
    
    if params.is_ensemble
        all_data.Uturb = all_Uturb;
        all_data.Vturb = all_Vturb;
        all_data.UturbVturb = all_UturbVturb;
        all_data.NanMask = all_NanMask;
    end
end

function weights = calculate_hanning_weights(cam_x, camera_no, overlap_info, camera_count)
    % Calculate Hanning weights for camera data based on overlap regions
    
    weights = ones(size(cam_x));
    
    if camera_no == 1 && camera_count > 1
        % First camera: apply right-side tapering
        if overlap_info.regions(camera_no) > 0
            overlap_size = overlap_info.regions(camera_no);
            window = 1 - overlap_info.hanning_windows{camera_no};
            weights(:, end-overlap_size+1:end) = weights(:, end-overlap_size+1:end) .* window;
        end
        
    elseif camera_no > 1 && camera_no < camera_count
        % Middle camera: apply both sides
        % Right overlap
        if overlap_info.regions(camera_no) > 0
            overlap_size = overlap_info.regions(camera_no);
            window = 1 - overlap_info.hanning_windows{camera_no};
            weights(:, end-overlap_size+1:end) = weights(:, end-overlap_size+1:end) .* window;
        end
        % Left overlap
        if overlap_info.regions(camera_no-1) > 0
            overlap_size = overlap_info.regions(camera_no-1);
            window = overlap_info.hanning_windows{camera_no-1};
            weights(:, 1:overlap_size) = weights(:, 1:overlap_size) .* window;
        end
        
    elseif camera_no > 1 && camera_no == camera_count
        % Last camera: apply left-side tapering
        if overlap_info.regions(camera_no-1) > 0
            overlap_size = overlap_info.regions(camera_no-1);
            window = overlap_info.hanning_windows{camera_no-1};
            weights(:, 1:overlap_size) = weights(:, 1:overlap_size) .* window;
        end
    end
end

function piv_result = interpolate_all_data(all_coords, all_data, coords, is_ensemble)
    % Single interpolation for all collected data
    
    % Remove NaN values
    valid_idx = ~isnan(all_data.u) & ~isnan(all_data.v) & ~isnan(all_coords.x) & ~isnan(all_coords.y);
    
    % Interpolate u velocity
    if sum(valid_idx) > 0
        interp = scatteredInterpolant(all_coords.x(valid_idx), all_coords.y(valid_idx), all_data.u(valid_idx), 'natural', 'none');
        piv_result.ux = interp(coords.x, coords.y);
        piv_result.ux(isnan(piv_result.ux)) = 0;
        
        % Interpolate v velocity
        interp = scatteredInterpolant(all_coords.x(valid_idx), all_coords.y(valid_idx), all_data.v(valid_idx), 'natural', 'none');
        piv_result.uy = interp(coords.x, coords.y);
        piv_result.uy(isnan(piv_result.uy)) = 0;
        
        % Interpolate b_mask
        interp = scatteredInterpolant(all_coords.x(valid_idx), all_coords.y(valid_idx), all_data.b_mask(valid_idx), 'natural', 'none');
        piv_result.b_mask = interp(coords.x, coords.y);
        piv_result.b_mask(isnan(piv_result.b_mask)) = 0;
        
        % Handle ensemble data
        if is_ensemble
            valid_turb = valid_idx & ~isnan(all_data.Uturb);
            if sum(valid_turb) > 0
                interp = scatteredInterpolant(all_coords.x(valid_turb), all_coords.y(valid_turb), all_data.Uturb(valid_turb), 'natural', 'none');
                piv_result.Uturb = interp(coords.x, coords.y);
                piv_result.Uturb(isnan(piv_result.Uturb)) = 0;
            else
                piv_result.Uturb = zeros(size(coords.x));
            end
            
            valid_turb = valid_idx & ~isnan(all_data.Vturb);
            if sum(valid_turb) > 0
                interp = scatteredInterpolant(all_coords.x(valid_turb), all_coords.y(valid_turb), all_data.Vturb(valid_turb), 'natural', 'none');
                piv_result.Vturb = interp(coords.x, coords.y);
                piv_result.Vturb(isnan(piv_result.Vturb)) = 0;
            else
                piv_result.Vturb = zeros(size(coords.x));
            end
            
            valid_turb = valid_idx & ~isnan(all_data.UturbVturb);
            if sum(valid_turb) > 0
                interp = scatteredInterpolant(all_coords.x(valid_turb), all_coords.y(valid_turb), all_data.UturbVturb(valid_turb), 'natural', 'none');
                piv_result.UturbVturb = interp(coords.x, coords.y);
                piv_result.UturbVturb(isnan(piv_result.UturbVturb)) = 0;
            else
                piv_result.UturbVturb = zeros(size(coords.x));
            end
            
            valid_mask = valid_idx & ~isnan(all_data.NanMask);
            if sum(valid_mask) > 0
                interp = scatteredInterpolant(all_coords.x(valid_mask), all_coords.y(valid_mask), all_data.NanMask(valid_mask), 'natural', 'none');
                piv_result.NanMask = interp(coords.x, coords.y);
                piv_result.NanMask(isnan(piv_result.NanMask)) = 0;
            else
                piv_result.NanMask = zeros(size(coords.x));
            end
        end
    else
        % No valid data
        piv_result.ux = zeros(size(coords.x));
        piv_result.uy = zeros(size(coords.x));
        piv_result.b_mask = zeros(size(coords.x));
        
        if is_ensemble
            piv_result.Uturb = zeros(size(coords.x));
            piv_result.Vturb = zeros(size(coords.x));
            piv_result.UturbVturb = zeros(size(coords.x));
            piv_result.NanMask = zeros(size(coords.x));
        end
    end
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
        % Build velocity file path with endpoint support
        if isempty(params.endpoint)
            vel_path = fullfile(params.calibrated_piv_dir, num2str(params.N), ...
                               ['Cam', num2str(camera_no)], params.data_subdir, ...
                               sprintf('%05d.mat', im_no));
        else
            vel_path = fullfile(params.calibrated_piv_dir, num2str(params.N), ...
                               ['Cam', num2str(camera_no)], params.data_subdir, params.endpoint, ...
                               sprintf('%05d.mat', im_no));
        end
        
        if ~exist(vel_path, 'file')
            error('Velocity file not found: %s', vel_path);
        end
        
        vel_data = load(vel_path);
        trunc = params.truncation(run_idx);
        
        % Extract and truncate data
        extracted_data = extract_camera_data(vel_data.piv_result(run_idx), trunc, params.is_ensemble);
        camera_data(camera_no).u = extracted_data.u;
        camera_data(camera_no).v = extracted_data.v;
        camera_data(camera_no).b_mask = extracted_data.b_mask;
        
        if params.is_ensemble
            camera_data(camera_no).Uturb = extracted_data.Uturb;
            camera_data(camera_no).Vturb = extracted_data.Vturb;
            camera_data(camera_no).UturbVturb = extracted_data.UturbVturb;
            camera_data(camera_no).NanMask = extracted_data.NanMask;
        end
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

function piv_result = post_process_masks(piv_result, is_ensemble)
    % Post-process binary masks
    
    piv_result.b_mask = piv_result.b_mask > 0.01;
    
    if is_ensemble
        piv_result.NanMask = piv_result.NanMask > 0.01;
    end
end





