function Gappy_POD_infil(setup,peakHeight,run,CameraNo)


directory_path = fullfile(setup.directory.base, 'CalibratedPIV', num2str(setup.imProperties.imageCount), ['Cam', num2str(CameraNo)], 'Instantaneous');
file_name = num2str(sprintf(setup.instantaneous.nameConvention{1}, 1));
file_path = fullfile(directory_path, file_name);
calibrated = load(file_path);

[rows, cols] = size(calibrated.piv_result(run).ux); % Grid size
U = zeros(rows* cols, setup.imProperties.imageCount); % Preallocate u
V = zeros(rows* cols, setup.imProperties.imageCount); % Preallocate v
masks = zeros(2*rows* cols, setup.imProperties.imageCount);
fprintf('Loading Gappy POD images at %s\n',(datetime('now')));
parfor imNo = 1:setup.imProperties.imageCount
    % Construct the filename
    file_name = num2str(sprintf(setup.instantaneous.nameConvention{1}, imNo));
    file_path = fullfile(directory_path, file_name);
    calibrated = load(file_path);
    unCalibrated = load(fullfile(setup.directory.base, 'UncalibratedPIV', num2str(setup.imProperties.imageCount), ['Cam', num2str(1)], 'Instantaneous', sprintf( setup.instantaneous.nameConvention{1},imNo)));
    % Store u and v in respective arrays (assuming rectangular grids)
    mask = (unCalibrated.piv_result(run).peak_choice ~=1) | (unCalibrated.piv_result(run).peak_mag < peakHeight);
    u= calibrated.piv_result(run).ux;
    u(mask) =0;
    v = calibrated.piv_result(run).uy;
    v(mask) =0;
    U(:,imNo) = u(:); % Store u component
    V(:,imNo) = v(:); % Store v component
    masks(:,imNo) = [mask(:);mask(:)];
end

% Reshape u and v into 2D matrices for decomposition
% Rows: spatial points (setup.imProperties.imageCount*m); Columns: time instances
masks = logical(masks);
X = [U;V];

for i = 1:25
    [U, S, V] = svd(X, 'econ'); % U: modes, S: singular values, V: temporal coefficients
    singularValues = diag(S);
    
    energy = singularValues.^2; % Energy is the square of the singular values
    totalEnergy = sum(energy); % Total energy
    
    cumulativeEnergy = cumsum(energy); 
    
    targetEnergy = 0.99 * totalEnergy;
    numModes = find(cumulativeEnergy >= targetEnergy, 1); % Find the first mode where the cumulative energy is >= 99%
    
    
    fprintf('The number of modes that take up 99%% of the energy is: %d\n', numModes);
 
    relevantModes = U(:,1:numModes);
    reconstruction_SVD = relevantModes * S(1:numModes, 1:numModes) * V(:, 1:numModes)';
    X(masks) = reconstruction_SVD(masks);
end


Lu = reshape(X(1:end/2, :), rows, cols, setup.imProperties.imageCount); % Low-rank u
Lv = reshape(X(end/2+1:end, :), rows, cols, setup.imProperties.imageCount); % Low-rank v
rejected_mask = reshape(masks(1:end/2,:),rows,cols,setup.imProperties.imageCount);




fprintf('Saving Gappy POD images at %s\n',(datetime('now')));
directory = fullfile(directory_path,'Gappy');
if ~exist(directory, 'dir')
    mkdir(directory);
end
parfor imNo = 1:setup.imProperties.imageCount
    % Construct the filename
    
    file_name = num2str(sprintf(setup.instantaneous.nameConvention{1}, imNo));
    file_path = fullfile(directory, file_name);
    unCalibrated = load(fullfile(setup.directory.base, 'UncalibratedPIV', num2str(setup.imProperties.imageCount), ['Cam', num2str(1)], 'Instantaneous', sprintf( setup.instantaneous.nameConvention{1},imNo)));
    b_mask = unCalibrated.piv_result(run).b_mask;

    piv_result=struct();
    piv_result(run).b_mask = b_mask;
    piv_result(run).ux = Lu(:,:,imNo);
    piv_result(run).uy = Lv(:,:,imNo);
    piv_result(run).Rejected_vectors = rejected_mask(:,:,imNo);
    PIVSAVE(file_path,piv_result)

    
end

fprintf('finished at %s\n',(datetime('now')));