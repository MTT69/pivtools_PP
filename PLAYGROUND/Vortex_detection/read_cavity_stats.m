       
base_dir = { ...
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

CameraNo = 1;
imageCount = 18000;

for base_idx = 1:length(base_dir)
    save_dir = fullfile(base_dir{base_idx}, 'Statistics', num2str(imageCount), ['Cam', num2str(CameraNo)], 'Instantaneous', 'Calibrated');
    save_path = fullfile(save_dir, 'cavity_stats.mat');
    load(save_path);
    % Loop through 11 images
    for img_idx = 1:11
        img_name = sprintf('Autocorrelation_convective velocity%d.jpg', img_idx);
        img_path = fullfile(save_dir, img_name);
        
        % Check if the file exists
        if exist(img_path, 'file')

            img = imread(img_path);
            pause(0.2)
            figure, imshow(img);
            title(sprintf('Image %d from %s', img_idx, base_dir{base_idx}), 'Interpreter', 'none');
        else
            warning('File not found: %s', img_path);
        end
    end
    close all
    convection_average = mean(convective_velocity(:,9))


    
end