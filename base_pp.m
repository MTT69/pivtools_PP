% Clear existing variables and close all figures
clear variables
close all

% Initialize the setup structure
setup = struct();

%%%%%%%%%%%% Global Processing Parameters %%%%%%%%%%%%%
% Define global type and endpoint for consistent processing
global_type = 'Instantaneous'; % Options: 'Instantaneous' or 'Ensemble' only
global_endpoint = ''; % Endpoint subfolder (e.g., 'rpca', 'denoised', etc.)
use_merged_data = false; % Set to true to process merged data instead of camera data

%%
%%%%%%%%%%%% Folders %%%%%%%%%%%%%
% Define setup directory structure
setup.directory = struct( ...
    'code', 'C:\Users\Lab8-2\Documents\pivtools_PP' ... % Location of PIV codes                
);

% base_dir = { ...
%    'D:\Full\Processed_PIV_validation\90degree_250light_250hz_1000dt'
% };

base_dir = { ...
   'D:\Full\Processed_PIV_validation\90degree_400light_100hz_3000dt'
};



run_instantaneous = [5];
run_ensemble = [5];
dt = {3000*10^(-6), 1000*10^(-6)}; % Time step in seconds for each run
frequency = {100,250};
setup.environment = struct( ...
    'local', true, ... % Specify execution environment: true for local, false for cluster (See compilation readme)     
    'numTasks', 3, ... % Number of parallel workers for heavy tasks                         
    'restartParpool', true, ... % Flag to reinitialize the parallel pool
    'imageLoadCores', 3, ... % Number of workers dedicated to image loading
    'maxCores', 3, ... % Maximum number of workers for light tasks
    'timeOut', 3600 ... % Timeout for idle workers in seconds
);



%%%%%%%%%%%% Setup File Loading %%%%%%%%%%%%%
% Configuration for loading existing setup files
setup.loadSetup = struct( ...
    'useMostRecent', true, ... % Use most recent setup file if true, otherwise use specific date
    'specificDate', '' ... % Specific date to load (format: yyyy_MM_dd)
);





%%%%%%%%%%%% Pipeline %%%%%%%%%%%%%
% Configure the processing pipeline settings
setup.pipeline = struct( ...
    'statistics_inst', false, ... % Generate instantaneous statistics
    'statistics_sum', false, ... % Generate instantaneous statistics
    'co_ords_gui', false, ... % edits coordinates interactively, set zero position
    'co_ords_complex', false, ... % edits coordinates, set bottom left value of each camera.
    'POST_POD', false,... % Perform POST POD analysis
    'pressure_reconstruction', false,... % Perform pressure reconstruction
    'manual_plots', false, ... % Generate manual plots
    'noise_floor', false, ... % Calculate noise floor
    'batch_BL_variation', false, ... % Perform batch boundary layer variation analysis
    'POD_rebuild', false, ... % Rebuild POD modes
    'vortex_video_gamma', false,...
    'gamma1_video', false, ... % Generate gamma1 video
    'gamma2_video', false, ... % Generate gamma1 video
    'lic_video', false,...
    'pressure_video', false,...
    'dot_probe', false,...
    'merge', false,...
    'flipUX',false,...
    'velocity_video', false,...
    'streamline_video', false...
);


%%%%%%%%% Figures %%%%%%%%%%%
% Define font sizes for figures
setup.figures.titleFontSize = 20; % Font size for titles
setup.figures.axisFontSize = 16; % Font size for axes
setup.figures.legendFontSize = 14; % Font size for legends
setup.figures.labelFontSize = 18; % Font size for labels



% Generate the base path including the parent directory
BasePath = genpath(setup.directory.code);
% Add the base path to MATLAB's search path
addpath(BasePath);
if isempty(gcp('nocreate'))
    Setup_parpool(setup, 'Images')
end
probe_x = 12;
probe_y = 0;


% plot_editor() can be used to edit any plots once opened.
CameraNo=1;
for i = 1:length(base_dir)
    setup.directory.base = base_dir{i};

    setup = load_setup_images(setup);
    setup.imProperties.dt =dt{i};
    setup.instantaneous.runs = run_instantaneous; % Number of runs to process
    setup.ensemble.runs = run_instantaneous; % Number of runs to process for ensemble statistics
    freq = frequency{i};
    if use_merged_data
        setup.imProperties.cameraCount = 1; % Merged data has only one camera
    end
    % Processing pipeline with global type and endpoint
    co_ord_editor(setup, CameraNo, global_type, global_endpoint, use_merged_data);
    multi_camera_coord_editor(setup, global_type, global_endpoint, use_merged_data)
    flip_ux_direction(setup, global_type, global_endpoint, use_merged_data)
    merge(setup, global_type, global_endpoint)
    Inst_statistics(setup, global_type, global_endpoint, use_merged_data)
    SumStatistics(setup, global_type, global_endpoint, use_merged_data)
    plot_maker_mean(setup, CameraNo, global_type, 'Test Data Analysis','xlabel', 'ylabel', global_endpoint, use_merged_data); % Create plots with custom title
    noise_floor(setup, CameraNo, global_type, global_endpoint, use_merged_data);
    POD_rebuild(setup, CameraNo, 'Below', global_type, use_merged_data);
    POST_POD_multi(setup, CameraNo, global_endpoint, 'Below', global_type, use_merged_data)
    POST_POD_multi(setup, CameraNo, global_endpoint, 'Above', global_type, use_merged_data)
    pressure_reconstruction(setup, CameraNo, global_type, global_endpoint, use_merged_data);
    BL_analysis(setup, CameraNo, global_endpoint, [-30,-5,30], global_type, use_merged_data)
    dot_probe(setup, CameraNo, global_endpoint, probe_x, probe_y, global_type, use_merged_data)
    Vortex_video_gamma(setup, global_endpoint, CameraNo, global_type, use_merged_data,freq)
    Gamma1_video(setup, global_endpoint, CameraNo, global_type, use_merged_data,freq)
    Gamma2_video(setup, global_endpoint, CameraNo, global_type, use_merged_data,freq)
    LIC_video(setup, global_endpoint, CameraNo, global_type, use_merged_data,freq)
    Pressure_video(setup, global_endpoint, CameraNo, global_type, use_merged_data,freq)
    make_velocity_video(setup, global_endpoint, CameraNo, global_type, use_merged_data, 'uy',freq)
    make_streamline_video(setup, global_endpoint, CameraNo, global_type, use_merged_data,freq)
    dot_probe_autocorr('D:\Full\Processed_PIV_validation\90degree_400light_100hz_3000dt\Statistics\18000\Cam1\Instantaneous\DotProbe\DotProbe_x12.00_y0.00_pass5.mat', freq)


end





