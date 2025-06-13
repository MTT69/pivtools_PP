function plot_save_mask(variable, mask,xcorners,ycorners,gcafontsize,titlefontsize,SaveLocation,runs_stats,wsize, loop, variableName,type,dir)
% plot_save_mask(variable, mask, xcorners, ycorners, gcafontsize, titlefontsize, SaveLocation, runs_stats, wsize, loop, variableName, type, dir)
%
% This function generates and saves a plot visualizing a given variable on a grid, with a corresponding mask 
% overlaid on top. The mask is used to hide certain areas of the variable, and a custom colormap is applied 
% based on the variable's type. The function supports different types of variables (e.g., 'Peak Heights', 'NaN Percentage', 
%  etc.) and automatically adjusts the display range and colormap based on the variable's properties.
%
% Inputs:
%   variable       - The data array (variable) to be visualized. This could represent various fields like velocity 
%                    or other physical properties over a grid.
%
%   mask           - A binary mask (same size as `variable`) used to hide certain areas in the plot by setting 
%                    those regions to `inf`. The mask will be overlaid on top of the variable visualization.
%
%   xcorners       - The x-coordinates for the grid to display the variable.
%
%   ycorners       - The y-coordinates for the grid to display the variable.
%
%   gcafontsize    - Font size for the axes and labels of the plot.
%
%   titlefontsize  - Font size for the title of the plot.
%
%   SaveLocation   - The directory where the figure will be saved as images in multiple formats (JPG, EPS, FIG).
%
%   runs_stats     - Array containing indices of the runs, used in the `limit` function to determine 
%                    appropriate bounds for the variable.
%
%   wsize          - Size of the window used for the PIV analysis, which will be used in the plot title.
%
%   loop           - The current loop index, used in conjunction with `runs_stats` to determine which data subset 
%                    is being processed.
%
%   variableName   - The name of the variable (e.g., 'Peak Heights', 'NaN Percentage').
%
%   type           - The type of data (e.g., 'Inst', 'Sum') used to determine whether the function should calculate 
%                    or load pre-saved displacement bounds.
%
%   dir            - Directory used for loading pre-saved displacement bounds if `type` is 'Sum'.
%
% Outputs:
%   This function does not return any outputs. It saves the generated plot in three formats (JPG, EPS, and FIG) to 
%   the specified `SaveLocation`.
%
% Description:
%   - The function first determines appropriate displacement limits (`lower_limit` and `upper_limit`) and colormap 
%     (`custommap`) based on the variable name (e.g., 'Peak Heights', 'NaN Percentage' etc.).
%     - For certain variable names like 'Peak Heights', 'NaN Percentage', or variables containing 'err', preset 
%       limits and colormaps are used.
%     - For other variables, the function calls the `limit` function to calculate the limits and determine the colormap.
%   - The variable array is modified by applying the mask (setting masked values to `inf`), and the data is plotted 
%     using `imagesc` with appropriate settings.
%   - A secondary `imagesc` plot is created for the mask, which is overlaid on top of the first plot. The second plot 
%     is hidden but ensures that the mask is correctly displayed in the final visualization.
%   - The plot is adjusted to maintain consistent aspect ratios and limits between the variable plot and the mask.
%   - The figure is saved in three formats (JPG, EPS, and FIG) to the `SaveLocation` directory. The file names include 
%     the `variableName` and window size (`wsize`).
%   - The figure is closed after saving to free up memory.
%
% Example:
%   % Plot and save the visualization for 'Assimilated U' with a custom mask, font sizes, and specified directory
%   plot_save_mask(velocity_data, mask, xcorners, ycorners, 12, 16, './plots', runs_stats, [32, 32], 1, 'Assimilated U', 'Inst', './');
%
% Notes:
%   - The mask is applied by creating a combined mask of original mask and NaN values, then setting those values to NaN.
%   - The `lower_limit` and `upper_limit` are chosen based on the type of variable (e.g., 'Assimilated U' or 'Peak Heights').
%   - The `imagesc` function is used for both the variable and the mask to create visual representations of the data. 
%     The `AlphaData` property ensures that only non-NaN values are displayed in the variable plot and masked areas
%     are shown in gray in the background layer.
%   - The figure is saved in three formats for compatibility with different software and use cases.
%   - The function is designed to handle various variable types by adjusting the colormap and bounds dynamically.
%   - LaTeX formatting is used for all text elements including tick labels, titles, and axis labels.
%


% Determine colormap and limits based on variable name
if contains(variableName, 'Peak Heights', 'IgnoreCase', true)
    lower_limit = 0;
    upper_limit = 1;
    custommap = 'parula';
elseif strcmp(variableName, 'NaN Percentage')
    lower_limit = 0;
    upper_limit = 100;
    custommap = 'parula';
elseif contains(variableName, 'err')
    lower_limit = -1;
    upper_limit = 1;
    custommap = redbluezero(lower_limit, upper_limit);
else
    % Use only valid (non-NaN and non-masked) data for limit calculation
    valid_data = variable(mask==0 & ~isnan(variable));
    [lower_limit, upper_limit, custommap] = limit(runs_stats, loop, valid_data, SaveLocation, variableName, type, dir);
    
    if lower_limit >= upper_limit || isinf(upper_limit) || lower_limit == upper_limit || (lower_limit <= 0 && upper_limit <= 0)
        lower_limit = -1;
        upper_limit = 1;
    end
end

% Create a combined mask: original mask OR NaN values
combined_mask = mask | isnan(variable);

% Apply combined mask by setting masked/NaN values to NaN
variable(combined_mask) = NaN;

% Convert combined mask to double for visualization
mask_vis = double(combined_mask);

% Create figure
figure('Visible', 'off')
fig = gcf;
set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

% Create variable data plot first (main layer)
ax1 = axes;
h = imagesc(xcorners, ycorners, variable, [lower_limit, upper_limit]);
set(h, 'AlphaData', ~isnan(variable)); % Ensures NaNs are fully transparent
set(gca, 'YDir', 'normal', 'FontSize', gcafontsize, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');



% Apply consistent colormap using othercolor('Spectral7') like plot_maker_instantaneous
if ischar(custommap)
    colormap(ax1, custommap);
else
    % Use the same colormap logic as plot_maker_instantaneous
    fullColormap = othercolor('Spectral7');
    blueToWhite = fullColormap(1:128,:);
    whiteToRed = fullColormap(129:end,:);
    
    if lower_limit >= 0
        % Only positive values: use white to red
        colormap(ax1, whiteToRed);
    elseif upper_limit <= 0
        % Only negative values: use blue to white
        colormap(ax1, blueToWhite);
    else
        % Both positive and negative: use full colormap
        colormap(ax1, fullColormap);
    end
end

% Hold on to add mask layer
hold(ax1, 'on');

% Create mask overlay with fixed light grey color (same as plot_maker_instantaneous)
[rows, cols] = size(mask_vis);
mask_rgb = cat(3, ones(rows, cols) * 0.8, ones(rows, cols) * 0.8, ones(rows, cols) * 0.8); % Light grey RGB

% Display the mask as an RGB image overlay
j = image(ax1, xcorners, ycorners, mask_rgb);
set(j, 'AlphaData', mask_vis * 0.7); % Show masked areas with transparency

colorbar;
daspect(ax1, [1 1 1]);
title(join([variableName, num2str(wsize(1,1)) 'x' num2str(wsize(1,2))], ""), 'FontSize', titlefontsize, 'Interpreter', 'latex');

hold(ax1, 'off');

% Save the figure in multiple formats
filename = fullfile(SaveLocation, join([variableName, num2str(wsize(1,1)), 'x', num2str(wsize(1,2))], ""));
saveas(fig, filename + ".jpg");
saveas(fig, filename + ".epsc");
saveas(fig, filename + ".fig");

% Close the figure and clean up
close(gcf);
delete(h);
delete(j);

end