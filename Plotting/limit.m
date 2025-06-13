
function [lower_limit, upper_limit, custommap]= limit(runs_stats, loop, variable,directory,variableName, type,dir)
    % [lower_limit, upper_limit, custommap] = limit(runs_stats, loop, variable, directory, variableName, type, dir)
    %
    % This function calculates the lower and upper displacement limits based on the provided data for graphing purposes, 
    % either using predefined bounds or calculating percentiles (1st and 99th). It saves these limits 
    % to a file for later use and generates a custom colormap for visualization based on the range of 
    % the variable values.
    %
    % Inputs:
    %   runs_stats    - Array containing statistics or indices for different runs. The first element is 
    %                   used to determine the first loop iteration.
    %
    %   loop          - The current loop index. This is used to determine whether it's the first loop 
    %                   (for saving displacement bounds) or subsequent loops (for loading saved bounds).
    %
    %   variable      - The data array for which the displacement bounds are to be calculated. This is 
    %                   used to compute the lower and upper percentiles.
    %
    %   directory     - The directory path where the displacement bounds file will be saved or loaded 
    %                   from.
    %
    %   variableName  - The name of the variable used for the displacement bounds file (e.g., 'ux', 'uy').
    %
    %   type          - The type of data being processed, which can be:
    %                   - 'Inst' or 'New' for instant or new data, where bounds are calculated and saved.
    %                   - 'Sum' for summed data, where bounds are loaded if previously saved.
    %                   - Any other type for loading previously saved bounds without recalculating.
    %
    %   dir           - Directory used to load the displacement bounds for the 'Sum' type.
    %
    % Outputs:
    %   lower_limit   - The calculated lower limit for the displacement based on the percentiles or loaded 
    %                   from a previously saved file.
    %
    %   upper_limit   - The calculated upper limit for the displacement based on the percentiles or loaded 
    %                   from a previously saved file.
    %
    %   custommap     - A custom colormap based on the range between the lower and upper displacement 
    %                   limits. The colormap is either `redbluezero` or `whitebluezero` depending on 
    %                   the values of `lower_limit`.
    %
    % Description:
    %   - For the first loop iteration (`loop == runs_stats(1)`), if the type is 'Inst' or 'New', the function 
    %     calculates the 1st and 99th percentiles of the variable and saves these as `lower_limit` and `upper_limit` 
    %     in a `.mat` file (`displacementbounds.mat`) in the specified directory.
    %   - If the type is 'Sum', the function attempts to load previously saved displacement bounds from the file 
    %     (`displacementbounds.mat`). If the file is not found, it falls back to calculating and saving new bounds
    %     this is so that the ensemble and instantaneous plots are on the same colourmap for easy comparison.
    %   - For any other type, the function loads the previously saved displacement bounds without recalculating them.
    %   - The function checks if `lower_limit` is less than 0 and, depending on this condition, generates a 
    %     custom colormap (`redbluezero` or `whitebluezero`).
    %
    % Example:
    %   % Calculate displacement bounds for 'ux' variable, save to directory, and generate a custom colormap
    %   [lower_limit, upper_limit, custommap] = limit(runs_stats, 1, ux_data, './results', 'ux', 'Inst', './');
    %
    % Notes:
    %   - The displacement bounds are saved in a `.mat` file under the specified `directory` with the name 
    %     `<variableName>displacementbounds.mat`.
    %   - The custom colormap (`custommap`) is created based on the `lower_limit` and `upper_limit` to visualize 
    %     the range of the displacement data.
    %   - The colormap is either `redbluezero` (if `lower_limit` is negative) or `whitebluezero` (if `lower_limit` 
    %     is non-negative).
    %   - The function uses the `prctile` function to calculate the 1st and 99th percentiles, ensuring the bounds 
    %     capture most of the data but exclude extreme outliers.
    %   - new constraint key word pairs can be added here easily for further control
    %


if loop == runs_stats(1) && strcmp(type, 'Inst') || strcmp(type, 'New')
    lower_limit = prctile(variable(:),1);
    upper_limit = prctile(variable(:),99);
    save(fullfile(directory, join([variableName,'displacementbounds.mat'],"")), 'lower_limit', 'upper_limit');
elseif strcmp(type, 'Sum')
    try
        load(fullfile(dir, join([variableName,'displacementbounds.mat'],"")));
    catch
        lower_limit = prctile(variable(:),1);
        upper_limit = prctile(variable(:),99);
        save(fullfile(directory, join([variableName,'displacementbounds.mat'],"")), 'lower_limit', 'upper_limit');
    end
else
    load(fullfile(directory, join([variableName,'displacementbounds.mat'],"")));
end

if lower_limit<0
    custommap = redbluezero(lower_limit,upper_limit);
    
else
    custommap = whitebluezero(lower_limit,upper_limit);
    
end

end