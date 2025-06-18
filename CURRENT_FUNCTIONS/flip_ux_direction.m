function flip_ux_direction(setup, data_type, endpoint, use_merged)
% FLIP_UX_DIRECTION - Multiplies all ux values by -1 in calibrated vector fields
%   This function processes all calibrated PIV data files and flips the sign
%   of the ux (x-velocity) component by multiplying by -1. The modified data
%   overwrites the original files.
%
%   Usage: flip_ux_direction(setup, CameraNo, data_type)
%
%   Inputs:
%       setup     - Setup structure containing file paths and configuration
%       CameraNo  - Camera number to process
%       data_type - 'Instantaneous' or 'Ensemble' data type
%
%   Example:
%       flip_ux_direction(setup, 1, 'Instantaneous');
%       flip_ux_direction(setup, 2, 'Ensemble');
%
%   WARNING: This function OVERWRITES the original files. Make sure to backup
%            your data before running this function!

    if nargin < 3
        endpoint = '';
    end
    if nargin < 4
        use_merged = false;
    end

    if setup.pipeline.flipUX
        fprintf('Flipping UX direction for type: %s at %s\n', data_type, datetime('now'));
        
        % Determine runs to process based on data type
        if strcmp(data_type, 'Instantaneous')
            runs_to_process = setup.instantaneous.runs;
            name_convention = setup.instantaneous.nameConvention{1};
            total_images = setup.imProperties.imageCount;
        elseif strcmp(data_type, 'Ensemble')
            runs_to_process = setup.ensemble.runs;
            name_convention = setup.ensemble.nameConvention{1};
            total_images = 1; % Ensemble typically has only one file per run
        else
            error('Invalid data_type. Must be "Instantaneous" or "Ensemble"');
        end
        
        for CameraNo = 1:setup.imProperties.cameraCount
            % Get data paths using centralized function
            paths = get_data_paths(setup, data_type, endpoint, CameraNo, use_merged);
            dataloc = paths.data_dir;
            
            if ~exist(dataloc, 'dir')
                warning('Data directory does not exist: %s', dataloc);
                continue;
            end
            

            % Process each run

            for run_idx = 1:length(runs_to_process)
                run = runs_to_process(run_idx);
                fprintf('Processing run %d/%d...\n', run_idx, length(runs_to_process));
                
                parfor img_no = 1:total_images 
                    % Construct filename
                    filename = sprintf(name_convention, img_no);
                    filepath = fullfile(dataloc, filename);
                
                    
                    try
                        % Load the PIV data
                        data = load(filepath);
                        
                        
                        % Get original ux values
                        original_ux = data.piv_result(run).ux;
                        
                        % Flip the sign of ux
                        data.piv_result(run).ux = -original_ux;
                        
                        
                        % Save the modified data back to the file
                        piv_result = data.piv_result;
                        PIVSAVE(filepath, piv_result);
                        
                        
                        % Progress indicator
                        if mod(img_no, 100) == 0 || img_no == total_images
                            fprintf('  Processed %d/%d files...\n', img_no, total_images);
                        end
                        
                    catch ME
                        warning('Error processing file %s: %s', filename, ME.message);
                        continue;
                    end
                end
                
                
            end
            fprintf('Completed flipping UX for camera %d\n', CameraNo);
        end
        
        fprintf('UX direction flipping completed at %s\n', datetime('now'));
    end
end
