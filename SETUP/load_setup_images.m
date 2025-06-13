function setup = load_setup_images(setup)

setupFiles = dir(fullfile(setup.directory.base, 'setup_*.mat'));

    if setup.loadSetup.useMostRecent
        % Sort by date and get most recent
        [~, idx] = max([setupFiles.datenum]);
        selectedFile = setupFiles(idx);
    else
        % Look for specific date
        targetPattern = sprintf('setup_%s_', setup.loadSetup.specificDate);
        matchIdx = contains({setupFiles.name}, targetPattern);
        if any(matchIdx)
            selectedFile = setupFiles(find(matchIdx, 1));
        else
            warning('Setup file for date %s not found in %s', setup.loadSetup.specificDate, setup.directory.base);
        end
    end
    
    % Load the setup file
    setup_load = load(fullfile(setup.directory.base, selectedFile.name));
    setup.imProperties = setup_load.setup.imProperties;
    setup.instantaneous = setup_load.setup.instantaneous;
    setup.ensemble = setup_load.setup.ensemble;


           

end