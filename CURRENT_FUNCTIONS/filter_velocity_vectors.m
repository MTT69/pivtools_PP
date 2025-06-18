function [u_filtered, v_filtered] = filter_velocity_vectors(u, v, mask, filter_type, varargin)
% FILTER_VELOCITY_VECTORS - Apply various filtering techniques to velocity vectors
%
% Inputs:
%   u, v - Velocity components
%   mask - Binary mask (true = invalid regions)
%   filter_type - 'gaussian', 'median', 'moving_average', 'outlier_removal'
%   varargin - Additional parameters depending on filter type
%
% Outputs:
%   u_filtered, v_filtered - Filtered velocity components

switch lower(filter_type)
    case 'gaussian'
        sigma = varargin{1}; % Gaussian standard deviation
        u_filtered = imgaussfilt(u, sigma, 'FilterDomain', 'spatial');
        v_filtered = imgaussfilt(v, sigma, 'FilterDomain', 'spatial');
        
    case 'median'
        kernel_size = varargin{1}; % e.g., [3 3]
        u_filtered = medfilt2(u, kernel_size);
        v_filtered = medfilt2(v, kernel_size);
        
    case 'moving_average'
        kernel_size = varargin{1}; % e.g., 3
        kernel = ones(kernel_size)/kernel_size^2;
        u_filtered = conv2(u, kernel, 'same');
        v_filtered = conv2(v, kernel, 'same');
        
    case 'outlier_removal'
        threshold = varargin{1}; % Standard deviation threshold
        
        % Calculate velocity magnitude
        vel_mag = sqrt(u.^2 + v.^2);
        
        % Remove outliers based on magnitude
        mean_mag = mean(vel_mag(~mask), 'omitnan');
        std_mag = std(vel_mag(~mask), 'omitnan');
        
        outlier_mask = vel_mag > (mean_mag + threshold * std_mag);
        
        u_filtered = u;
        v_filtered = v;
        u_filtered(outlier_mask) = NaN;
        v_filtered(outlier_mask) = NaN;
        
        % Apply median filter to fill gaps
        u_filtered = fillmissing(u_filtered, 'movmedian', 5);
        v_filtered = fillmissing(v_filtered, 'movmedian', 5);
        
    otherwise
        error('Unknown filter type: %s', filter_type);
end

% Apply mask
u_filtered(mask) = NaN;
v_filtered(mask) = NaN;

end
