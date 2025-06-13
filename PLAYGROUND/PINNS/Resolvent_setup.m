base = { ...
    'D:\Full\Processed_PIV_validation\90degree_400light_100hz_3000dt', ...
%     'D:\Full\Processed_PIV_validation\90degree_250light_250hz_1000dt', ...
%     'D:\Full\Processed_PIV_validation\60degree_400light_100hz_3000dt_reverse', ...
%     'D:\Full\Processed_PIV_validation\60degree_250light_250hz_1000dt_reverse', ...
%     'D:\Full\Processed_PIV_validation\60degree_400light_100hz_3000dt', ...
%     'D:\Full\Processed_PIV_validation\60degree_250light_250hz_1000dt', ...
%     'D:\Full\Processed_PIV_validation\30degree_400light_100hz_3000dt_reverse', ...
%     'D:\Full\Processed_PIV_validation\30degree_250light_250hz_1000dt_reverse', ...
%     'D:\Full\Processed_PIV_validation\30degree_400light_100hz_3000dt', ...
%     'D:\Full\Processed_PIV_validation\30degree_250light_250hz_1000dt', ...
};

imageCount = 18000;
caseImages = 3600;
CameraNo =1;
endpoint = "";
run = 5;
aperature = 17;
speeds = [0.35, 1.18,0.35, 1.18,0.35, 1.18,0.35, 1.18,0.35, 1.18,];

for i = 1:length(base)
    statsloc = fullfile(base{i}, 'Statistics', num2str(imageCount), ['Cam', num2str(CameraNo)], 'Instantaneous', 'Calibrated');
    % Load Co_ords.mat
    VelData = load(fullfile(statsloc, 'MeanStats16x16.mat'));
    mean_U = VelData.mean_U;
    mean_V = VelData.mean_V;
    dataloc = (fullfile(base{i}, 'CalibratedPIV', num2str(imageCount),['Cam', num2str(CameraNo)], 'Instantaneous',endpoint));
    VelData = load(fullfile(dataloc, [num2str(sprintf('%05d.mat', 1))]));
    Co_ords = load(fullfile(dataloc, 'Co_ords.mat'));
    b_mask = VelData.piv_result(run).b_mask;
    uncalibrated_dataloc = (fullfile(base{i}, 'UncalibratedPIV', num2str(imageCount),['Cam', num2str(CameraNo)], 'Instantaneous','Corr'));
    correlation_stats = load(fullfile(uncalibrated_dataloc,'5.mat'));

    PINN = fullfile(base{i}, 'PINN', num2str(imageCount), ['Cam', num2str(CameraNo)], 'Instantaneous', 'Calibrated');

    if ~exist(PINN, 'dir')
        mkdir(PINN);
    end

    datamask = (correlation_stats.CorMean < 0.7) | b_mask;
    Re = aperature*10^-3

% 
%     need_reynolds number
%     need co-ords
    



end