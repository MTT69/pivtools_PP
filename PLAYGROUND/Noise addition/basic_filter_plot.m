% filepath: /c:/Users/Lab8-2/Documents/PIVTOOLS/PIVTOOLS/PLAYGROUND/Noise addition/basic_filter_plot.m

% Define the root folder
rootFolder = 'D:\Synthetic_processed\Images_5snr_15_mean_2px';

% Define paths to each type of filtered image
originalPath = fullfile(rootFolder, 'BasicFilters', 'Cam1', 'Batch1', 'original');
ssminPath = fullfile(rootFolder, 'BasicFilters', 'Cam1', 'Batch1', 'ssmin');
timePath = fullfile(rootFolder, 'BasicFilters', 'Cam1', 'Batch1', 'time');

% Load frame 1 images
original_frame1 = imread(fullfile(originalPath, 'filtered-original''-frame1.tiff'));
ssmin_frame1 = imread(fullfile(ssminPath, 'filtered-ssmin-frame1.tiff'));
time_frame1 = imread(fullfile(timePath, 'filtered-time-frame1.tiff'));

% Load frame 2 images
original_frame2 = imread(fullfile(originalPath, 'filtered-original''-frame2.tiff'));
ssmin_frame2 = imread(fullfile(ssminPath, 'filtered-ssmin-frame2.tiff'));
time_frame2 = imread(fullfile(timePath, 'filtered-time-frame2.tiff'));

% Convert to double for consistent processing
original_frame1 = double(original_frame1);
ssmin_frame1 = double(ssmin_frame1);
time_frame1 = double(time_frame1);
original_frame2 = double(original_frame2);
ssmin_frame2 = double(ssmin_frame2);
time_frame2 = double(time_frame2);

% Extract central 128x128 pixels from each image
[rows, cols] = size(original_frame1);
centerRow = round(rows/2);
centerCol = round(cols/2);
halfSize = 64; % Half of 128

rowRange = (centerRow-halfSize):(centerRow+halfSize-1);
colRange = (centerCol-halfSize):(centerCol+halfSize-1);

original_frame1_center = original_frame1(rowRange, colRange);
ssmin_frame1_center = ssmin_frame1(rowRange, colRange);
time_frame1_center = time_frame1(rowRange, colRange);
original_frame2_center = original_frame2(rowRange, colRange);
ssmin_frame2_center = ssmin_frame2(rowRange, colRange);
time_frame2_center = time_frame2(rowRange, colRange);

% Find common colorbar limits for all images
minVal = min([
    original_frame1_center(:); ssmin_frame1_center(:); time_frame1_center(:);
    original_frame2_center(:); ssmin_frame2_center(:); time_frame2_center(:)
]);
maxVal = max([
    original_frame1_center(:); ssmin_frame1_center(:); time_frame1_center(:);
    original_frame2_center(:); ssmin_frame2_center(:); time_frame2_center(:)
]);

% Create figure for Frame 1
figure('Position', [100 100 1200 300]);

% Frame 1 subplot
subplot(1, 3, 1);
imagesc(original_frame1_center, [minVal maxVal]);
title('Original - Frame 1');
axis image;
colorbar;
colormap gray;

subplot(1, 3, 2);
imagesc(time_frame1_center, [minVal maxVal]);
title('Time - Frame 1');
axis image;
colorbar;
colormap gray;

subplot(1, 3, 3);
imagesc(ssmin_frame1_center, [minVal maxVal]);
title('SSMin - Frame 1');
axis image;
colorbar;
colormap gray;

% Save Frame 1 figure
saveas(gcf, fullfile(rootFolder, 'BasicFilters', 'Cam1', 'Batch1', 'frame1_plot.png'));

% Create figure for Frame 2
figure('Position', [100 500 1200 300]);

% Frame 2 subplot
subplot(1, 3, 1);
imagesc(original_frame2_center, [minVal maxVal]);
title('Original - Frame 2');
axis image;
colorbar;
colormap gray;

subplot(1, 3, 2);
imagesc(time_frame2_center, [minVal maxVal]);
title('Time - Frame 2');
axis image;
colorbar;
colormap gray;

subplot(1, 3, 3);
imagesc(ssmin_frame2_center, [minVal maxVal]);
title('SSMin - Frame 2');
axis image;
colorbar;
colormap gray;

% Save Frame 2 figure
saveas(gcf, fullfile(rootFolder, 'BasicFilters', 'Cam1', 'Batch1', 'frame2_plot.png'));