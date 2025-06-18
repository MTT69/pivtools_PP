function plot_editor()
% PLOT_EDITOR Interactive editor for loaded .fig files
%   This function provides an interface to edit figure properties with
%   configurable text/font settings in code and interactive colorbar bounds.
%   
%   Usage:
%   1. Load a .fig file: openfig('your_figure.fig')
%   2. Edit the configuration variables below
%   3. Run: plot_editor()
%   4. Adjust colorbar bounds interactively
%
%   Features:
%   - Configurable labels and font sizes in code
%   - Interactive colorbar bounds adjustment
%   - Real-time preview of changes
%   - Option to save edited figure

%% CONFIGURATION SECTION - EDIT THESE VALUES
% Set update_labels to true to apply the label changes below
update_labels = true;
figure(1)
% Label settings (set to empty string '' to keep current value)
new_title = 'My Custom Title';
new_xlabel = 'X Position (mm)';
new_ylabel = 'Y Position (mm)';
new_cblabel = 'Velocity (m/s)';

% Font size settings (set to 0 to keep current value)
new_title_fontsize = 16;
new_axis_fontsize = 12;
new_label_fontsize = 14;

%% MAIN FUNCTION - DO NOT EDIT BELOW THIS LINE
% Check if there's a current figure
if isempty(gcf) || ~isvalid(gcf)
    error('No figure found. Please load a .fig file first using openfig(''filename.fig'')');
end

fig = gcf;
ax = gca;

fprintf('=== PLOT EDITOR ===\n');
fprintf('Current figure: Figure %d\n\n', fig.Number);

% Get current values for comparison and fallback
current_title = get(get(ax, 'Title'), 'String');
current_xlabel = get(get(ax, 'XLabel'), 'String');
current_ylabel = get(get(ax, 'YLabel'), 'String');
current_title_fontsize = get(get(ax, 'Title'), 'FontSize');
current_axis_fontsize = get(ax, 'FontSize');
current_label_fontsize = get(get(ax, 'XLabel'), 'FontSize');

% Get colorbar if it exists
cb = findobj(fig, 'Type', 'colorbar');
current_cblabel = '';
if ~isempty(cb)
    current_cblabel = get(get(cb, 'Label'), 'String');
end

% Apply label and font changes if requested
if update_labels
    fprintf('APPLYING LABEL AND FONT UPDATES:\n');
    
    % Use current values if new ones are empty/zero
    if isempty(new_title)
        final_title = current_title;
    else
        final_title = new_title;
    end
    
    if isempty(new_xlabel)
        final_xlabel = current_xlabel;
    else
        final_xlabel = new_xlabel;
    end
    
    if isempty(new_ylabel)
        final_ylabel = current_ylabel;
    else
        final_ylabel = new_ylabel;
    end
    
    if isempty(new_cblabel)
        final_cblabel = current_cblabel;
    else
        final_cblabel = new_cblabel;
    end
    
    if new_title_fontsize == 0
        final_title_fontsize = current_title_fontsize;
    else
        final_title_fontsize = new_title_fontsize;
    end
    
    if new_axis_fontsize == 0
        final_axis_fontsize = current_axis_fontsize;
    else
        final_axis_fontsize = new_axis_fontsize;
    end
    
    if new_label_fontsize == 0
        final_label_fontsize = current_label_fontsize;
    else
        final_label_fontsize = new_label_fontsize;
    end
    
    % Apply changes
    title(ax, final_title, 'FontSize', final_title_fontsize, 'Interpreter', 'latex');
    xlabel(ax, final_xlabel, 'FontSize', final_label_fontsize, 'Interpreter', 'latex');
    ylabel(ax, final_ylabel, 'FontSize', final_label_fontsize, 'Interpreter', 'latex');
    set(ax, 'FontSize', final_axis_fontsize);
    
    if ~isempty(cb) && ~isempty(final_cblabel)
        ylabel(cb, final_cblabel, 'FontSize', final_label_fontsize, 'Interpreter', 'latex');
    end
    
    fprintf('✓ Title: %s (Font: %.1f)\n', char(final_title), final_title_fontsize);
    fprintf('✓ X-label: %s (Font: %.1f)\n', char(final_xlabel), final_label_fontsize);
    fprintf('✓ Y-label: %s (Font: %.1f)\n', char(final_ylabel), final_label_fontsize);
    fprintf('✓ Axis font size: %.1f\n', final_axis_fontsize);
    if ~isempty(cb)
        fprintf('✓ Colorbar label: %s (Font: %.1f)\n', char(final_cblabel), final_label_fontsize);
    end
else
    fprintf('Label updates disabled (update_labels = false)\n');
end

% Colorbar bounds adjustment
fprintf('\nCOLORBAR BOUNDS ADJUSTMENT:\n');
if ~isempty(cb)
    current_clim = get(ax, 'CLim');
    fprintf('Current colorbar bounds: [%.6f, %.6f]\n', current_clim(1), current_clim(2));
    
    fprintf('Enter new bounds (press Enter with empty string to finish):\n');
    
    while true
        % Get lower bound
        lower_input = input(sprintf('Lower bound (current: %.6f): ', current_clim(1)), 's');
        if isempty(lower_input)
            fprintf('Colorbar editing finished.\n');
            break;
        end
        
        lower_bound = str2double(lower_input);
        if isnan(lower_bound)
            fprintf('Invalid input. Please enter a number.\n');
            continue;
        end
        
        % Get upper bound
        upper_input = input(sprintf('Upper bound (current: %.6f): ', current_clim(2)), 's');
        if isempty(upper_input)
            fprintf('Colorbar editing finished.\n');
            break;
        end
        
        upper_bound = str2double(upper_input);
        if isnan(upper_bound)
            fprintf('Invalid input. Please enter a number.\n');
            continue;
        end
        
        % Validate bounds
        if lower_bound >= upper_bound
            fprintf('Lower bound must be less than upper bound. Try again.\n');
            continue;
        end
        
        % Apply new bounds
        set(ax, 'CLim', [lower_bound, upper_bound]);
        current_clim = [lower_bound, upper_bound];
        
        fprintf('Applied bounds: [%.6f, %.6f]\n', lower_bound, upper_bound);
        fprintf('Enter new bounds or press Enter with empty string to finish:\n');
    end
else
    fprintf('No colorbar found in the current figure.\n');
end

% Option to save the edited figure
fprintf('\nSAVE OPTIONS:\n');
save_fig = false;
while true
    response = input('Save the edited figure? (y/n): ', 's');
    if strcmpi(response, 'y') || strcmpi(response, 'yes')
        save_fig = true;
        break;
    elseif strcmpi(response, 'n') || strcmpi(response, 'no') || isempty(response)
        save_fig = false;
        break;
    else
        fprintf('Please enter y or n\n');
    end
end

if save_fig
    % Get the current figure's filename if it exists
    figName = get(fig, 'FileName');
    
    if ~isempty(figName)
        % Extract filename without extension
        [filepath, filename, ~] = fileparts(figName);
        default_name = filename;
        default_path = filepath;
    else
        % No existing filename, use current directory
        default_name = 'edited_figure';
        default_path = pwd;
    end
    
    % Ask for filename with current name as default
    filename = input(sprintf('Enter filename (default: %s): ', default_name), 's');
    if isempty(filename)
        filename = default_name;
    end
    
    % Remove extension if provided
    [~, filename, ~] = fileparts(filename);
    
    % Create full path
    if ~isempty(default_path)
        full_filename = fullfile(default_path, filename);
    else
        full_filename = filename;
    end
    
    % Save in multiple formats, overwriting existing files
    try
        saveas(fig, [full_filename, '.fig']);
        saveas(fig, [full_filename, '.jpg']);
        saveas(fig, [full_filename, '.epsc']);
        fprintf('Figure overwritten/saved as:\n');
        fprintf('  %s.fig\n', full_filename);
        fprintf('  %s.jpg\n', full_filename);
        fprintf('  %s.epsc\n', full_filename);
    catch ME
        fprintf('Error saving figure: %s\n', ME.message);
    end
end

fprintf('\nPlot editing completed!\n');

end
