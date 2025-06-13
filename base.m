% Clear existing variables and close all figures
clear variables
close all

% Initialize the setup structure
setup = struct();
%%
%%%%%%%%%%%% Folders %%%%%%%%%%%%%
% Define setup directory structure
setup.directory = struct( ...
    'code', 'C:\Users\Lab8-2\Documents\PIVTOOLS\PIVTOOLS', ... % Location of PIV codes       
    'python', '', ... % Location of Python executable for IM7/IMS reading                       
    'calibration', '', ... % Location of pinhole calibration file (if required)                 
    'gslCluster', '/local/software/gsl/2.7_gcc13', ... % Location of GSL source files           
    'fftwCluster', '/local/software/fftw/3.3.10' ... % Location of FFTW source files            
);

%%
% Define other required paths
base_dir = { ...
   'D:\Synthetic_processed\Images_20snr_15_mean_2px'
};

source = { ...
    'D:\Synthetic_raw\images_20_snr_15_bias'
};  % Cell array of source file roots; images should reside within a "Cam1" folder        
% NotD: IMS source should be the path to the source .set file, not just the folder. ExamplD: '/data/PIVCam1.set'
% NotD: if using multiple cases, e.g save the source suffixes as _0, _1 but leave source with no suffix this is applied automatically

uinf = [0.35, 1.18,0.35, 1.18,0.35, 1.18,0.35, 1.18,0.35, 1.18,];
sample = [100,250,100,250,100,250,100,250,100,250]; %hz 
dt = {0.0075}; % basic calibration dt - cell for different runs


%%
%%%%%%%%%%% Environment %%%%%%%%%%%%%


setup.environment = struct( ...
    'local', true, ... % Specify execution environment: true for local, false for cluster      
    'numTasks', 6, ... % Number of parallel workers for heavy tasks                         
    'restartParpool', true, ... % Flag to reinitialize the parallel pool
    'imageLoadCores', 6, ... % Number of workers dedicated to image loading
    'maxCores', 6 ... % Maximum number of workers for light tasks
);


%%%%%%%%%%%% Image Properties %%%%%%%%%%%%%
setup.imProperties = struct( ...
    'caseImages', 1000, ... % Number of image pairs in one run of PIV                       
    'imageCount',1000, ... % Total number of images (including repeats)                
    'batchSize', 1000, ... % Number of image pairs processed in a single batch               
    'parforbatch', 50, ... % Number of images to process in a sub-pass                      
    'imageSize', [1025, 1025], ... % Image dimensions: [height, width]                      
    'imageType', 'im', ... % Image format options: 'im', 'im7', 'ims', 'cine'
    'timeResolved', false, ... % Flag for time-resolved imaging
    'cameraCount', 1, ... % Number of cameras used
    'combineRuns', false, ... % Process repeats as one if true (looks for _0, _1 folders)
    'scaleFactor',3.416666667, ... % Calibration scale factor (pixels per mm)
    'yOffset', 0, ... % Calibration offset for the coordinate system
    'xOffset',  0 ...
);



% xOffset = {1234,1256, 1246}; % a normal b reversed c 90
% degree
% xOffset = {1246}; % basical calibration x offset 


%% Instantaneous naming conventions

% Naming conventions for images to be loaded in
setup.imProperties.nameConvention = {...
    'B%05d_A.tif', ...                     % Convention 1: 'B#####_a.tiff' %%% ensure the suffix matches the file type
    'B%05d_B.tif'};

% Naming conventions for processedImages to be Saved
setup.imProperties.saveNameConvention = {...
'B%05d_a.tiff', ...  % Convention 1: 'B#####_a.tiff'
'B%05d_b.tiff'};     % Convention 2: 'B#####_b.tiff'
%%
%%%%%%%%%%%% Pipeline %%%%%%%%%%%%%
% Configure the processing pipeline settings
setup.pipeline = struct( ...
    'createMask', false, ... % Create a mask for processing
    'loadMask', false, ... % Load a pre-existing mask
    'polygonsToRemove', 1, ... % Polygons to remove or use [top, right, bottom, left] for a square mask
    'instantaneous', true, ... % Run instantaneous PIV
    'ensemble', true, ... % Run ensemble averaging
    'calibrate_inst', false, ... % Enable calibration
    'calibrate_sum', true, ... % Enable calibration
    'calibrateType', 'basic', ... % Type of calibration method ('basic', etc.)
    'compile', false, ... % Compile C files (forced true on cluster)
    'medianFilter', false, ... % Apply median filtering for outlier detection
    'globalFilter', false, ... % Apply global filtering (standard deviation approach)
    'statistics_inst', false, ... % Generate instantaneous statistics
    'statistics_correlation', false, ... % Generate correlation statistics
    'statistics_sum', true, ... % Generate instantaneous statistics
    'prefilter', false, ... % Save pre-filtered images outside PIV workstream
    'storePlanes', true, ... % Store ensemble correlation planes
    'calculateSumWindow', false, ... % Automatically calculate the optimal sum window size
    'RPCA', false,... % Perform RPCA infilling
    'RPCA_split', false,... % Perform RPCA infilling on top and bottom of image
    "Fluctuation_videos",false,...
    'SPOD', false,... % Perform SPOD analysis
    'POST_POD', false,... % Perform POST POD analysis
    'POST_POD_stats', false,... % Perform POST POD analysis
    'POST_POD_batches', false...
);



%%
%%%%%%%%%%%% Instantaneous PIV %%%%%%%%%%%%%
% Configure parameters for instantaneous PIV
% Configure parameters for instantaneous PIV
setup.instantaneous = struct( ...
    'windowSize', [128 128; 64 64; 32 32; 16 16; 16 16; 16 16], ... % Correlation window sizes
    'overlap', [50, 50, 50, 50, 50, 50], ... % Percent overlap for each window level
    'dt', 1, ... % Time steps between images to load (e.g., 2 loads images 1 and 3 as a pair)
    'peakFinder', 'gauss6', ... % Method for peak finding ('gauss3', etc.)
    'windowType', 'gaussian', ... % Type of correlation window ('blackman', 'gaussian')
    'runs', [6], ... % List of runs used for PIV analysis
    'interpolator', 'spline', ... % Interpolation method ('spline', 'linear')
    'secondaryPeaks', false... % use secondary or tertiary peaks if first peak in an outlier
);
setup.instantaneous.nameConvention = {'%05d.mat'}; % Naming convention for saving processed vectors


%%
%%%%%%%%%%%% Ensemble PIV %%%%%%%%%%%%%
% Configure parameters for ensemble PIV
setup.ensemble = struct( ...
    'windowSize', [64 64; 32 32; 16 16; 4 4], ... % Correlation window sizes
    'resumeCase', 0,... % resume ensemble from previous run % 0 means start fresh
    'overlap', [0, 50, 50, 50], ... % Overlap percentages for each window level
    'sumWindow', [16, 16], ... % Bounding window size for ensemble, dynamically settable
    'runs', [3,4], ... % List of runs used for PIV analysis
    'interpolator', 'spline', ... % Interpolation method ('spline', 'linear')
    'windowType', 'square', ... % Window shape ('square'), as described in literature
    'dt', 1, ... % Time steps between images to load (e.g., 2 loads images 1 and 3 as a pair)
    'convergedRun', 4 ... % Run number for calculating SumWindow for a single pixel
);
setup.ensemble.type = {'std', 'std', 'std', 'single', 'single'}; % Averaging types ('std' or 'single')
setup.ensemble.nameConvention  = {'%05d.mat'}; % Naming convention for saving ensemble results

%%
%%%%%%%%%%%% Filters %%%%%%%%%%%%%
% filters can be found in loading_and_pre_processing with explanations for setup
% Time and POD filters recommended
filters = cell(1, 1); 

ssmin_filt_size = [3 3];
gauss_filt_size = [3 3];
norm_filt_size = [3 3];
median_filt_size = [3 3];
filters{4} = struct('type', 'gaussian', 'size', gauss_filt_size, 'sigma', 0.5);
filters{1} = struct('type', 'time', 'size', 50);
filters{2} = struct('type', 'ssmin', 'size', ssmin_filt_size);
% filters{1} = struct('type', 'norm', 'size', norm_filt_size, 'max_gain', 1);
filters{3} = struct('type', 'median', 'size', median_filt_size);
%%
%%%%%%%%% Figures %%%%%%%%%%%
% Define font sizes for figures
setup.figures.titleFontSize = 20; % Font size for titles
setup.figures.axisFontSize = 16; % Font size for axes
setup.figures.legendFontSize = 14; % Font size for legends
setup.figures.labelFontSize = 18; % Font size for labels


%%


% Generate the base path including the parent directory
BasePath = genpath(setup.directory.code);
% Add the base path to MATLAB's search path
addpath(BasePath);
% Automatically set up additional parameters
[setup, python_environment] = Auto_setup(setup);

CameraNo = 1;

for i = 1:length(source)
    setup.directory.base = base_dir{i};
    if ~exist(setup.directory.base, 'dir')
       mkdir(setup.directory.base)
    end
    u = uinf(i);
    SR = sample(i);
    save(fullfile(setup.directory.base, sprintf('setup_%s_%d_images.mat', datetime('now', 'Format', 'yyyy, MM, dd'), setup.imProperties.imageCount)), 'setup');
    setup.directory.source = source{i};
    setup.imProperties.dt =dt{i};
    masks = Masking(setup); % sets up mask object
    [setup, filters]= Perform_pre_filter(setup, filters,masks);  % TODO to leave this here filters needs to be set per casethis filters and saves all images before running anything and sets the source directory to the filtered images
   
    Perform_PIV(setup, masks, filters); % perform instantaneous PIV
    Perform_PIV_calibration_inst(setup, CameraNo,'Instantaneous',setup.instantaneous.runs,setup.imProperties.imageCount)
    Correlation_stats(setup) % calculates efficacy of instantaneous PIV
    Inst_statistics(setup,'Calibrated',CameraNo,'')


%     POST_POD_UV(setup,'Calibrated',CameraNo,'','Below')
    POST_POD_UV(setup,'Calibrated',CameraNo,'','Above')
    POST_POD_UV_Batches(setup,'Calibrated',CameraNo,'','Below')
%     POST_POD_UV_Batches(setup,'Calibrated',CameraNo,'','Above')
    Frequency_post_pod(setup,'Calibrated',CameraNo,'','Below', 100, 20, [1,2,3,4,5],0.3)
    Fluctuation_videos(setup, '', u, SR, 'v', 'Below', [1,2,3,4,5], 75, 20, CameraNo, 'Calibrated')%Fluctuation_videos(setup, endpoint, u_inf, sampling_rate, plot_type, domain, modes, dividing_row, frame_range)

    RPCA_infil(setup, 0.02, 100,1e-7, 0,true,CameraNo)                    %(setup, lambda, maxIter,tol, peakHeight,gappy,CameraNo)
    RPCA_infil_top_bottom(setup, 0.02, 100,1e-7, 0,true,CameraNo,80) %RPCA_infil_top_bottom(setup, lambda, maxIter,tol, peakHeight, gappy, CameraNo, split_row)


    SPOD_stats(setup,5,CameraNo,'',512,true)                                % SPOD_stats(setup,run,CameraNo,endpoint,nDFT,combine)
    reconstructSPOD_modes(setup,CameraNo,'',true, 1:20, 2:3)                % reconstructSPOD_modes(setup,CameraNo,endpoint,combine, modes, frequencies)
    reconstruct_temporal_SPOD(setup,CameraNo,'',true, 1:20, 2:3,512,5)      % reconstruct_temporal_SPOD(setup,CameraNo,endpoint,combine, modes, frequencies,nDFT,run)
    plot_SPOD_reconstruction(setup, CameraNo, '', true, 1:20, 2:3, 500,5)   % plot_SPOD_reconstruction(setup, CameraNo, endpoint, combine, modes, frequencies, numTimeElements,run)



    %%%%%%%%%%%%%%%%%Ensemble PIV%%%%%%%%%%%%%%%%%%%%%%%
    
    EnsemblePIV(setup, masks,filters) % perform sum of correlation method PIV
    Perform_PIV_calibration_sum(setup, CameraNo,'Ensemble',setup.ensemble.runs)
    SumStatistics(setup,'Calibrated',CameraNo)
end





