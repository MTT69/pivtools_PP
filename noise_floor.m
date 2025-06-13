function noise_floor(setup, CameraNo)
    for i = setup.instantaneous.runs
        directory=fullfile(setup.directory.base, 'Statistics', num2str(setup.imProperties.imageCount), ['Cam' num2str(CameraNo)], 'Instantaneous','Calibrated');

        dataloc = (fullfile(setup.directory.base, 'CalibratedPIV', num2str(setup.imProperties.imageCount),['Cam', num2str(CameraNo)], 'Instantaneous'));
        MeanStats=load(fullfile(directory,['MeanStats',num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2)),'.mat']));


        mean_u_prime_squared_calibrated = MeanStats.U_prime_Uprime_mean;
        mean_v_prime_squared_calibrated = MeanStats.V_prime_Vprime_mean;
 
        Co_ords = load(fullfile(dataloc, 'Co_ords.mat'));
        VelData = load(fullfile(dataloc, [num2str(sprintf(setup.instantaneous.nameConvention{1}, 1))]));
        b_mask = VelData.piv_result(i).b_mask;
        ycorners = [Co_ords.Co_ords(i).y(1,1), Co_ords.Co_ords(i).y(end, end)]; %, xcoords(m, 1), xcoords(m, n)];
        xcorners = [Co_ords.Co_ords(i).x(1,1), Co_ords.Co_ords(i).x(end, end)];
        
        
        % Calculate pixel-wise noise field
        uncertainty = 0.1/(setup.imProperties.scaleFactor/10^(-3)*setup.imProperties.dt);
        noise_field = ones(size(mean_u_prime_squared_calibrated)) * uncertainty;
        
        % Calculate noise percentage for each pixel
        total_signal = sqrt(mean_u_prime_squared_calibrated + mean_v_prime_squared_calibrated);
        noise_percentage_field = 100 * (noise_field ./ (total_signal));
        
        % Plot noise percentage field
        variableName = 'Noise Percentage';
        plot_save_mask(noise_percentage_field, b_mask, xcorners, ycorners, setup.figures.axisFontSize, setup.figures.titleFontSize, directory, setup.instantaneous.runs, [setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)], i, variableName, 'Inst', 'dir');

    end
end