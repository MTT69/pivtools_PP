% Define the data location and parameters
clear variables
dataloc = 'D:\Full\Processed_PIV_validation\90degree_400light_100hz_3000dt\UncalibratedPIV\18000\Cam1\Instantaneous';
batchSize = 3600;      % number of images per batch
totalImages = 18000;
numBatches = totalImages / batchSize;  % should be 5 batches

i = 5;                 % index to select the appropriate PIV_result element

% Preallocate matrix for the averaged profiles (each column is 77x1 for one batch)
profiles = zeros(77, numBatches);

% Loop over batches
for b = 1:numBatches
    % Initialize an accumulator for the profile over this batch
    profile_sum = zeros(77, 1);
    
    % Loop over images within the current batch
    for imNo = (b-1)*batchSize + 1 : b*batchSize
        disp(imNo)
        % Construct file name with leading zeros (e.g., '00001.mat')
        fileName = sprintf('%05d.mat', imNo);
        % Load the velocity data file
        VelData = load(fullfile(dataloc, fileName));
        % Extract the U velocity component from the specified PIV result index
        U = VelData.piv_result(i).ux;
        
        % Extract the region of interest:
        %   - rows 1:77 correspond to y direction
        %   - columns 10:15 correspond to the x vector between the 10th and 15th position
        region = U(1:77, 10:15);
        
        % Average the extracted region along the x-direction (columns) to form a 77x1 profile
        image_profile = mean(region, 2);
        
        % Accumulate the profile for this batch
        profile_sum = profile_sum + image_profile;
    end
    
    % Compute the average profile for the current batch
    profiles(:, b) = profile_sum / batchSize;
end

% Plotting: Swap axes, improve aesthetics, and use LaTeX labels
figure;
hold on;
colors = lines(numBatches); % Generate distinguishable colors

for b = 1:numBatches
    plot((1:77)', profiles(:, b), 'Color', colors(b, :), 'LineWidth', 2, ...
         'DisplayName', sprintf('Batch %d', b));
end

% LaTeX-formatted labels
xlabel('\textbf{y Positions (Vectors Above Wall)}', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('\textbf{U Displacement (Pixels)}', 'Interpreter', 'latex', 'FontSize', 14);
title('\textbf{Batch Profile Comparison}', 'Interpreter', 'latex', 'FontSize', 16);

% Improve legend visibility
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);

% Enhance plot aesthetics
grid on;
box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5); % Adjust axes appearance
hold off;

