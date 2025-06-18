% Clear existing variables and close all figures
clear variables
close all

% Initialize the setup structure
setup = struct();
%%
%%%%%%%%%%%% Folders %%%%%%%%%%%%%
% Define setup directory structure
setup.directory = struct( ...
    'code', 'C:\Users\Lab8-2\Documents\pivtools_PP' ... % Location of PIV codes                
);

base_dir = { ...
   'D:\Full\Processed_PIV_validation\90degree_250light_250hz_1000dt'
};

run_instantaneous = [5];
run_ensemble = [6];
dt = {1000*10^(-6), 3000*10^(-6)}; % Time step in seconds for each run

setup.environment = struct( ...
    'local', true, ... % Specify execution environment: true for local, false for cluster (See compilation readme)     
    'numTasks', 6, ... % Number of parallel workers for heavy tasks                         
    'restartParpool', true, ... % Flag to reinitialize the parallel pool
    'imageLoadCores', 6, ... % Number of workers dedicated to image loading
    'maxCores', 6, ... % Maximum number of workers for light tasks
    'timeOut', 3600 ... % Timeout for idle workers in seconds
);



%%%%%%%%%%%% Setup File Loading %%%%%%%%%%%%%
% Configuration for loading existing setup files
setup.loadSetup = struct( ...
    'useMostRecent', true, ... % Use most recent setup file if true, otherwise use specific date
    'specificDate', '2024_01_15' ... % Specific date to load (format: yyyy_MM_dd)
);





%%%%%%%%%%%% Pipeline %%%%%%%%%%%%%
% Configure the processing pipeline settings
setup.pipeline = struct( ...
    'statistics_inst', false, ... % Generate instantaneous statistics
    'statistics_sum', false, ... % Generate instantaneous statistics
    'ensemble_cords', false, ... % edits ensemble coordinates
    'instantaneous_cords', false, ... % edits instantaneous coordinates
    'POST_POD', false,... % Perform POST POD analysis
    'pressure_reconstruction', false,... % Perform pressure reconstruction
    'manual_plots', false, ... % Generate manual plots
    'noise_floor', false, ... % Calculate noise floor
    'batch_BL_variation', true, ... % Perform batch boundary layer variation analysis
    'POD_rebuild', false, ... % Rebuild POD modes
    'vortex_video_gamma', false,...
    'gamma1_video', false, ... % Generate gamma1 video
    'gamma2_video', false, ... % Generate gamma1 video
    'lic_video', false,...
    'pressure_video', false,...
    'dot_probe', true...
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

CameraNo = 1;

% plot_editor() can be used to edit any plots once opened.

for i = 1:length(base_dir)
    setup.directory.base = base_dir{i};

    setup = load_setup_images(setup);
    setup.imProperties.dt =dt{i};
    setup.instantaneous.runs = run_instantaneous; % Number of runs to process
    setup.ensemble.runs = run_instantaneous; % Number of runs to process for ensemble statistics
    co_ord_editor(setup, CameraNo); % Edit coordinates for the camera
    plot_maker_instantaneous(setup, CameraNo, 'Calibrated', 'Test Data Analysis','xlabel', 'ylabel'); % Create plots with custom title
    Inst_statistics(setup,'Calibrated',CameraNo,'')
    SumStatistics(setup,'Calibrated',CameraNo)
    noise_floor(setup, CameraNo) % Calculate noise floor for the camera
    POD_rebuild(setup, CameraNo, false); % Rebuild POD modes for the camera
    pressure_reconstruction(setup, CameraNo);
    POST_POD_multi(setup, CameraNo, '', 'Below')
    POST_POD_multi(setup, CameraNo, '', 'Above')
    BL_analysis(setup, CameraNo, '', [-30,-5,30])
    Vortex_video_gamma(setup, '' ,CameraNo)
    Gamma1_video(setup, '', CameraNo)
    Gamma2_video(setup, '', CameraNo)
    LIC_video(setup, '', CameraNo)
    Pressure_video(setup, '', CameraNo)
    dot_probe(setup, CameraNo, '', probe_x, probe_y)







    % 
  % plot_SPOD_reconstruction(setup, CameraNo, endpoint, combine, modes, frequencies, numTimeElements,run)





end





