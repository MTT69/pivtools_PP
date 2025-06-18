function paths = get_data_paths(setup, type, endpoint, CameraNo, use_merged)
    % GET_DATA_PATHS - Centralized function to determine data and statistics paths
    % Inputs:
    %   setup - setup structure
    %   type - 'Instantaneous' or 'Ensemble'
    %   endpoint - endpoint subfolder (e.g., 'rpca', 'denoised', etc.) or ''
    %   CameraNo - camera number (ignored for merged data)
    %   use_merged - logical, true to use merged data, false for traditional camera data
    % Outputs:
    %   paths - structure containing data_dir and stats_dir
    
    if nargin < 5
        use_merged = false; % Default to traditional camera data
    end
    
    base_dir = setup.directory.base;
    
    if use_merged
        % Merged data structure: Merged/type/endpoint/
        if isempty(endpoint)
            paths.data_dir = fullfile(base_dir, 'Merged', type);
            paths.stats_dir = fullfile(base_dir, 'Statistics', 'Merged', type);
        else
            paths.data_dir = fullfile(base_dir, 'Merged', type, endpoint);
            paths.stats_dir = fullfile(base_dir, 'Statistics', 'Merged', type, endpoint);
        end
        paths.pressure_dir = fullfile(paths.data_dir, 'Pressure');
        paths.video_dir = fullfile(base_dir, 'Videos', 'Merged', type);
        
    else
        % Traditional camera data structure: CalibratedPIV/imageCount/CamX/type/endpoint/
        if isempty(endpoint)
            paths.data_dir = fullfile(base_dir, 'CalibratedPIV', num2str(setup.imProperties.imageCount), ...
                                     ['Cam', num2str(CameraNo)], type);
        else
            paths.data_dir = fullfile(base_dir, 'CalibratedPIV', num2str(setup.imProperties.imageCount), ...
                                     ['Cam', num2str(CameraNo)], type, endpoint);
        end
        
        % Statistics always go to Calibrated subdirectory for traditional data
        paths.stats_dir = fullfile(base_dir, 'Statistics', num2str(setup.imProperties.imageCount), ...
                                  ['Cam', num2str(CameraNo)], type, 'Calibrated');
        paths.pressure_dir = fullfile(paths.data_dir, 'Pressure');
        paths.video_dir = fullfile(base_dir, 'Videos', num2str(setup.imProperties.imageCount), ...
                                  ['Cam', num2str(CameraNo)]);
    end
end
