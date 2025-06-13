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
   'D:\Full\Processed_PIV_validation\90degree_400light_100hz_3000dt'
};

run_instantaneous = [5];
run_ensemble = [6];
dt = {3000*10^(-6), 1000*10^(-6)}; % Time step in seconds for each run

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
    'RPCA', false,... % Perform RPCA infilling
    'RPCA_split', false,... % Perform RPCA infilling on top and bottom of image
    'SPOD', false,... % Perform SPOD analysis
    'POST_POD', false,... % Perform POST POD analysis
    'pressure_reconstruction', false... % Perform pressure reconstruction
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



for i = 1:length(base_dir)
    setup.directory.base = base_dir{i};

    setup = load_setup_images(setup);
    setup.imProperties.dt =dt{i};
    setup.instantaneous.runs = run_instantaneous; % Number of runs to process
    setup.ensemble.runs = run_instantaneous; % Number of runs to process for ensemble statistics
%     co_ord_editor(setup, CameraNo); % Edit coordinates for the camera
    % plot_maker_instantaneous(setup, CameraNo, 'Calibrated', 'Test Data Analysis','xlabel', 'ylabel'); % Create plots with custom title
    Inst_statistics(setup,'Calibrated',CameraNo,'')
%     noise_floor(setup, CameraNo) % Calculate noise floor for the camera
    % Cavity-only reconstruction (default)
    % POD_rebuild(setup, 1);
    

    
    pressure_reconstruction(setup, CameraNo);



    % SumStatistics(setup,'Calibrated',CameraNo)


    % POST_POD_UV(setup,'Calibrated',CameraNo,'','Above')
    % POST_POD_UV_Batches(setup,'Calibrated',CameraNo,'','Below')
    % Frequency_post_pod(setup,'Calibrated',CameraNo,'','Below', 100, 20, [1,2,3,4,5],0.3)
    % Fluctuation_videos(setup, '', u, SR, 'v', 'Below', [1,2,3,4,5], 75, 20, CameraNo, 'Calibrated')%Fluctuation_videos(setup, endpoint, u_inf, sampling_rate, plot_type, domain, modes, dividing_row, frame_range)

    % RPCA_infil(setup, 0.02, 100,1e-7, 0,true,CameraNo)                    %(setup, lambda, maxIter,tol, peakHeight,gappy,CameraNo)
    % RPCA_infil_top_bottom(setup, 0.02, 100,1e-7, 0,true,CameraNo,80) %RPCA_infil_top_bottom(setup, lambda, maxIter,tol, peakHeight, gappy, CameraNo, split_row)


    % SPOD_stats(setup,5,CameraNo,'',512,true)                                % SPOD_stats(setup,run,CameraNo,endpoint,nDFT,combine)
    % reconstructSPOD_modes(setup,CameraNo,'',true, 1:20, 2:3)                % reconstructSPOD_modes(setup,CameraNo,endpoint,combine, modes, frequencies)
    % reconstruct_temporal_SPOD(setup,CameraNo,'',true, 1:20, 2:3,512,5)      % reconstruct_temporal_SPOD(setup,CameraNo,endpoint,combine, modes, frequencies,nDFT,run)
    % plot_SPOD_reconstruction(setup, CameraNo, '', true, 1:20, 2:3, 500,5)   % plot_SPOD_reconstruction(setup, CameraNo, endpoint, combine, modes, frequencies, numTimeElements,run)





end





