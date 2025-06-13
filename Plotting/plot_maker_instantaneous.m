function plot_maker_instantaneous(setup, CameraNo, Type, customTitle, xLabel, yLabel)
% PLOT_MAKER_INSTANTANEOUS Creates instantaneous plots from PIV data with mask functionality
%   Inputs:
%       setup       - setup configuration structure
%       cameraNo    - camera number identifier
%       Type        - 'Calibrated' or 'Uncalibrated' data type
%       customTitle - optional custom title prefix for plots
%       xLabel      - optional custom x-axis label (default: 'X Position')
%       yLabel      - optional custom y-axis label (default: 'Y Position')

    % Set default Type if not provided
    if nargin < 3
        Type = 'Calibrated';
    end
    
    % Set default custom title if not provided
    if nargin < 4 || isempty(customTitle)
        customTitle = '';
    end
    
    % Set default axis labels if not provided
    if nargin < 5 || isempty(xLabel)
        xLabel = 'X Position';
    end
    
    if nargin < 6 || isempty(yLabel)
        yLabel = 'Y Position';
    end
    
    % Set endpoint (empty for now, can be extended later)
    endpoint = '';
    
    % Determine data location based on Type
    if strcmp(Type,'Calibrated')        
        dataloc = (fullfile(setup.directory.base, 'CalibratedPIV', num2str(setup.imProperties.imageCount),['Cam', num2str(CameraNo)], 'Instantaneous',endpoint));
    elseif strcmp(Type, 'Uncalibrated')
        dataloc = (fullfile(setup.directory.base, 'UncalibratedPIV', num2str(setup.imProperties.imageCount),['Cam', num2str(CameraNo)], 'Instantaneous'));
    else
        error('Incompatible type for instantaneous statistics');
    end
    
    directory=fullfile(setup.directory.base, 'Statistics', num2str(setup.imProperties.imageCount), ['Cam' num2str(CameraNo)], 'Instantaneous',Type);
    
    selectedField = '';  % Store user selection
    
    for i = setup.instantaneous.runs
        % Load velocity data and coordinates
        VelData = load(fullfile(dataloc, [num2str(sprintf(setup.instantaneous.nameConvention{1}, 1))]));
        Co_ords = load(fullfile(dataloc, 'Co_ords.mat'));
        
        % Extract mask from velocity data
        b_mask = VelData.piv_result(i).b_mask;
        
        MeanStats=load(fullfile(directory,['MeanStats',num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2)),'.mat']));
        
        % On first iteration, show available plots and get user selection
        if i == setup.instantaneous.runs(1)
            fieldNames = fieldnames(MeanStats);
            fprintf('Available plots:\n');
            for j = 1:length(fieldNames)
                fprintf('%d - %s\n', j, fieldNames{j});
            end
            
            while true
                userChoice = input('Select plot to display (enter number): ');
                if userChoice >= 1 && userChoice <= length(fieldNames) && mod(userChoice,1) == 0
                    selectedField = fieldNames{userChoice};
                    fprintf('Selected: %s\n', selectedField);
                    break;
                else
                    fprintf('Invalid selection. Please enter a number between 1 and %d.\n', length(fieldNames));
                end
            end
            
            % Get user input for bounds
            fprintf('\nEnter custom bounds (press Enter to skip for auto-calculation):\n');
            lowerBoundInput = input('Lower bound: ', 's');
            upperBoundInput = input('Upper bound: ', 's');
            
            if ~isempty(lowerBoundInput) && ~isempty(upperBoundInput)
                userLowerBound = str2double(lowerBoundInput);
                userUpperBound = str2double(upperBoundInput);
                if ~isnan(userLowerBound) && ~isnan(userUpperBound)
                    useCustomBounds = true;
                    fprintf('Using custom bounds: [%.3f, %.3f]\n', userLowerBound, userUpperBound);
                else
                    useCustomBounds = false;
                    fprintf('Invalid bounds entered. Using auto-calculation.\n');
                end
            else
                useCustomBounds = false;
                fprintf('Using auto-calculation for bounds.\n');
            end
        end
        
        ycorners = [Co_ords.Co_ords(i).y(1,1), Co_ords.Co_ords(i).y(end, end)];
        xcorners = [Co_ords.Co_ords(i).x(1,1), Co_ords.Co_ords(i).x(end, end)];
        
        % Get the selected field data
        plotData = MeanStats.(selectedField);
        
        % Apply mask functionality using plot_save_mask approach
        % Create a combined mask: original mask OR NaN values
        combined_mask = b_mask | isnan(plotData);
        
        % Apply combined mask by setting masked/NaN values to NaN
        plotData(combined_mask) = NaN;
        
        % Convert combined mask to double for visualization
        mask_vis = double(combined_mask);
        
        % Create figure with mask functionality
        figure('Name', sprintf('%s - Window %dx%d', selectedField, setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)));
        
        % Determine bounds using valid data only
        if useCustomBounds
            lowerLimit = userLowerBound;
            upperLimit = userUpperBound;
        else
            % Use only valid (non-NaN and non-masked) data for limit calculation
            valid_data = plotData(~combined_mask & ~isnan(plotData));
            
            if isempty(valid_data)
                lowerLimit = -1;
                upperLimit = 1;
            else
                % Smart colormap based on data range
                dataMin = min(valid_data(:));
                dataMax = max(valid_data(:));
                
                if dataMin >= 0
                    lowerLimit = 0;
                    upperLimit = prctile(valid_data(:), 99);
                elseif dataMax <= 0
                    lowerLimit = prctile(valid_data(:), 1);
                    upperLimit = 0;
                else
                    maxAbsValue = max(abs([prctile(valid_data(:), 1), prctile(valid_data(:), 99)]));
                    lowerLimit = -maxAbsValue;
                    upperLimit = maxAbsValue;
                end
            end
        end
        
        % Create variable data plot first (main layer)
        ax1 = axes;
        h = imagesc(xcorners, ycorners, plotData, [lowerLimit, upperLimit]);
        set(h, 'AlphaData', ~isnan(plotData)); % Ensures NaNs are fully transparent
        set(gca, 'YDir', 'normal', 'FontSize', setup.figures.axisFontSize, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');
        
        % Store original data for data cursor functionality
        setappdata(h, 'OriginalData', MeanStats.(selectedField));
        setappdata(h, 'Mask', combined_mask);
        
        % Apply colormap based on data characteristics using othercolor('Spectral7')
        fullColormap = othercolor('Spectral7');
        blueToWhite = fullColormap(1:128,:);
        whiteToRed = fullColormap(129:end,:);
        
        if lowerLimit >= 0
            % Only positive values: use white to red
            colormap(ax1, whiteToRed);
        elseif upperLimit <= 0
            % Only negative values: use blue to white
            colormap(ax1, blueToWhite);
        else
            % Both positive and negative: use full colormap
            colormap(ax1, fullColormap);
        end
        
        % Hold on to add mask layer
        hold(ax1, 'on');
        
        % Create mask overlay with fixed light grey color
        % Create a light grey RGB image for masked areas
        [rows, cols] = size(mask_vis);
        mask_rgb = cat(3, ones(rows, cols) * 0.8, ones(rows, cols) * 0.8, ones(rows, cols) * 0.8); % Light grey RGB
        
        % Display the mask as an RGB image overlay
        j = image(ax1, xcorners, ycorners, mask_rgb);
        set(j, 'AlphaData', mask_vis * 0.7); % Show masked areas with transparency
        
        % Set up custom data cursor to show original velocity values
        dcm = datacursormode(gcf);
        set(dcm, 'UpdateFcn', @(obj, event_obj) customDataCursor(obj, event_obj, h));
        
        colorbar;
        daspect(ax1, [1 1 1]);
        
        % Create title with optional custom prefix
        if ~isempty(customTitle)
            titleText = customTitle;
        else
            titleText = sprintf('%s - Window Size: %dx%d', selectedField, setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2));
        end
        
        title(titleText, 'FontSize', setup.figures.titleFontSize, 'Interpreter', 'latex');
        xlabel(xLabel, 'FontSize', setup.figures.labelFontSize, 'Interpreter', 'latex');
        ylabel(yLabel, 'FontSize', setup.figures.labelFontSize, 'Interpreter', 'latex');
        
        % Ensure consistent figure aspect ratio
        daspect([1 1 1]);
        
        hold(ax1, 'off');
    end
end

function txt = customDataCursor(~, event_obj, imageHandle)
    % Custom data cursor function to show original velocity values
    
    % Get position
    pos = get(event_obj, 'Position');
    
    % Get original data and mask from stored app data
    originalData = getappdata(imageHandle, 'OriginalData');
    mask = getappdata(imageHandle, 'Mask');
    
    % Get image data properties
    imageData = get(imageHandle, 'CData');
    xData = get(imageHandle, 'XData');
    yData = get(imageHandle, 'YData');
    
    % Calculate pixel indices
    if length(xData) == 2
        % XData and YData are limits
        xIdx = round((pos(1) - xData(1)) / (xData(2) - xData(1)) * (size(imageData, 2) - 1)) + 1;
        yIdx = round((pos(2) - yData(1)) / (yData(2) - yData(1)) * (size(imageData, 1) - 1)) + 1;
    else
        % XData and YData are full coordinate arrays
        [~, xIdx] = min(abs(xData - pos(1)));
        [~, yIdx] = min(abs(yData - pos(2)));
    end
    
    % Ensure indices are within bounds
    xIdx = max(1, min(size(originalData, 2), xIdx));
    yIdx = max(1, min(size(originalData, 1), yIdx));
    
    % Get values
    originalValue = originalData(yIdx, xIdx);
    isMasked = mask(yIdx, xIdx);
    
    % Create display text
    if isMasked
        txt = {['X: ', num2str(pos(1), '%.3f')], ...
               ['Y: ', num2str(pos(2), '%.3f')], ...
               'Value: MASKED'};
    else
        txt = {['X: ', num2str(pos(1), '%.3f')], ...
               ['Y: ', num2str(pos(2), '%.3f')], ...
               ['Value: ', num2str(originalValue, '%.6f')]};
    end
end