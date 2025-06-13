base = { ...
    'D:\Full\Processed_PIV_validation\90degree_400light_100hz_3000dt', ...
    'D:\Full\Processed_PIV_validation\90degree_250light_250hz_1000dt', ...
    'D:\Full\Processed_PIV_validation\60degree_400light_100hz_3000dt_reverse', ...
    'D:\Full\Processed_PIV_validation\60degree_250light_250hz_1000dt_reverse', ...
    'D:\Full\Processed_PIV_validation\60degree_400light_100hz_3000dt', ...
    'D:\Full\Processed_PIV_validation\60degree_250light_250hz_1000dt', ...
    'D:\Full\Processed_PIV_validation\30degree_400light_100hz_3000dt_reverse', ...
    'D:\Full\Processed_PIV_validation\30degree_250light_250hz_1000dt_reverse', ...
    'D:\Full\Processed_PIV_validation\30degree_400light_100hz_3000dt', ...
    'D:\Full\Processed_PIV_validation\30degree_250light_250hz_1000dt', ...
};
samplingRate =[100,250,100,250,100,250,100,250,100,250,]; %hz 
dt = {3000*10^(-6),}; % basic calibration dt - cell for different runs
min_expected_freq =0.3;
imageCount = 18000;
caseImages = 3600;
endpoint = '';
run = 5;
x_loc = -5;
scaleFactor = 9.53;
CameraNo =1;
close all


for i = 1:length(base)
        dataloc = fullfile(base{i}, 'Statistics', num2str(imageCount), ['Cam', num2str(CameraNo)], 'Instantaneous', 'Calibrated');
        % Load Co_ords.mat
        VelData = load(fullfile(dataloc, 'MeanStats16x16.mat'));
        mean_ux_calibrated = VelData.mean_U;
        mean_u_prime_squared_calibrated = VelData.U_prime_Uprime_mean;
        mean_v_prime_squared_calibrated = VelData.V_prime_Vprime_mean;

        mean_ux = mean_ux_calibrated *((scaleFactor/10.^-3*dt{i}));
        mean_u_prime_squared = mean_u_prime_squared_calibrated *((scaleFactor/10.^-3*dt{i}))^2;
        mean_u_prime = sqrt(mean_u_prime_squared);

        mean_v_prime_squared = mean_v_prime_squared_calibrated * ((scaleFactor / 10.^-3 * dt{i}))^2;
        mean_v_prime = sqrt(mean_v_prime_squared);

        
        figure(1)
        imagesc(mean_u_prime, [0,1]);
        ylim([75,130])
        xlim([100,225])
        daspect([1 1 1])
        colorbar
        title('Average Absolute U Prime',FontSize=16)
        
        figure(2)
        imagesc(mean_v_prime, [0,1]);
        ylim([75,130])
        xlim([100,225])
        daspect([1 1 1])
        colorbar
        title('Average Absolute V Prime',FontSize=16)



        close(gcf)
end