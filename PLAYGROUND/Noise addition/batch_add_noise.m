% Script to batch process images with gaussian noise
% This script adds noise to all images and saves them in a new folder

% Parameters for noise
snr = 50;       % SNR value
bias = 15;     % Mean noise value

% Set up paths
inputFolder = 'D:\Synthetic_raw\images_0_noise\Cam1';
outputFolder = 'D:\Synthetic_raw\images_50_snr_15_bias\Cam1';

% Create output folder if it doesn't exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
    fprintf('Created output directory: %s\n', outputFolder);
end

% Get list of all image files
fileList = dir(fullfile(inputFolder, 'B*_*.tif'));
fileNames = {fileList.name};

% Process each file
totalFiles = numel(fileNames);
fprintf('Found %d files to process\n', totalFiles);

% Start progress display
tic;
fprintf('Processing images: 0%%');

parfor i = 1:totalFiles
    % Load image
    fileName = fileNames{i};
    inputPath = fullfile(inputFolder, fileName);
    I = imread(inputPath);
    
    % Apply noise
    noisyI = gaussian_noise(I, snr, bias);
    noisyI = uint8(noisyI);
    
    % Save processed image
    outputPath = fullfile(outputFolder, fileName);
    imwrite(noisyI, outputPath);
    
    % Display progress
    if mod(i, 50) == 0 || i == totalFiles
        fprintf('\rProcessing images: %.1f%%', i/totalFiles*100);
    end
end

% Finish progress display
elapsedTime = toc;
fprintf('\nProcessing completed in %.2f seconds\n', elapsedTime);
fprintf('All images saved to: %s\n', outputFolder);
