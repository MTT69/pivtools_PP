function Setup_parpool(setup, type)
% Setup_parpool(setup, type)
%
% Set up a parallel pool for processing tasks using either process-based or 
% thread-based parallel execution. The function checks if a parallel pool 
% already exists and, if not, creates one based on the specified type.
%
% Inputs:
%   setup - Struct containing environment setup parameters, including:
%           - environment.numTasks: Number of workers (processes or threads) to use.
%           - environment.timeOut: Idle timeout setting for the parallel pool.
%
%   type  - String specifying the type of parallel execution. Options are:
%           - 'Processes': Parallel pool using multiple processes.
%           - 'Threads': Parallel pool using multiple threads.
%
% Description:
%   - The function first checks if a parallel pool is already running using 
%     `gcp('nocreate')`. If a pool exists, it does nothing.
%   - If no pool is running, it sets up a new parallel pool based on the 
%     specified `type`. For 'Processes', it launches a pool with a specified 
%     number of processes, while for 'Threads', it launches a thread-based pool.
%   - If the type is 'Processes', the `IdleTimeout` property is set according 
%     to the `setup.environment.timeOut` parameter, which controls how long the 
%     pool remains open when idle.
%   - If an invalid type is provided, the function throws an error and provides 
%     guidance on valid options.
%
% Example:
%   % Set up a process-based parallel pool with the specified number of workers
%   Setup_parpool(setup, 'Processes');
%
% Notes:
%   - `parpool('local', ...)` is used for both process-based and thread-based 
%     environments, but the underlying execution model differs depending on the 
%     type specified.
%   - Ensure that the `numTasks` field in the `setup` struct is set appropriately 
%     based on the available CPU cores and the task requirements.
%   - The `IdleTimeout` setting ensures that the parallel pool does not remain 
%     active indefinitely when idle.


    % Set up parallel pool for processing
    

    
    % Configure the parallel pool based on the specified type
    switch type
        case 'Processes'
            % Process-based parallel execution
            c = parcluster('local');
            c.NumWorkers = setup.environment.numTasks;
            saveProfile(c);
            numProcesses = setup.environment.numTasks; % Number of processes
            
            % Check if a parallel pool already exists
            poolObj = gcp('nocreate');  % Get current pool, but don't create new one
        
            % If no pool exists or the number of workers doesn't match the desired count, create a new pool
            if isempty(poolObj) || poolObj.NumWorkers ~= numProcesses
                delete(poolObj);  % Delete existing pool if the worker count is wrong
                parpool('local', numProcesses, 'IdleTimeout', setup.environment.timeOut);  % Create a new pool
                fprintf('Created a new process-based pool with %d workers.\n', numProcesses);
            else
                fprintf('Using existing process-based pool with %d workers.\n', poolObj.NumWorkers);
            end
        case 'Images'
            c = parcluster('local');
            c.NumWorkers = setup.environment.imageLoadCores; % Specify number of workers
            saveProfile(c);
            numProcesses = setup.environment.imageLoadCores; % Number of processes
            
            % Check if a parallel pool exists and if the number of workers matches
            poolObj = gcp('nocreate');  % Get the current pool without creating a new one
            
            % Only create a new pool if it doesn't exist or has a different number of workers
            if isempty(poolObj) || poolObj.NumWorkers ~= numProcesses
                delete(poolObj);  % Delete the existing pool if it doesn't match the desired size
                parpool('local', numProcesses, 'IdleTimeout', setup.environment.timeOut);  % Create a new pool
                fprintf('Created a new process-based environment for image loading with %d workers.\n', numProcesses);
            else
                fprintf('Using the existing process-based environment with %d workers.\n', poolObj.NumWorkers);
            end
        case 'Max'
            c = parcluster('local');
            c.NumWorkers = setup.environment.maxCores; % Specify number of workers
            saveProfile(c);
            numProcesses = setup.environment.maxCores; % Number of processes
            
            % Check if a parallel pool exists and if the number of workers matches
            poolObj = gcp('nocreate');  % Get the current pool without creating a new one
            
            % Only create a new pool if it doesn't exist or has a different number of workers
            if isempty(poolObj) || poolObj.NumWorkers ~= numProcesses
                delete(poolObj);  % Delete the existing pool if it doesn't match the desired size
                parpool('local', numProcesses, 'IdleTimeout', setup.environment.timeOut);  % Create a new pool
                fprintf('Created a new process-based environment for max processing with %d workers.\n', numProcesses);
            else
                fprintf('Using the existing process-based environment with %d workers.\n', poolObj.NumWorkers);
            end


            
        
            
        otherwise
            % Error handling for invalid types
            error('Invalid parallel setup specified');
    end
end