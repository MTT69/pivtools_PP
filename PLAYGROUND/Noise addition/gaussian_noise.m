% Example usage with SNR (higher SNR = less noise)
function noisyImage = gaussian_noise(I, snr, bias)
    % GAUSSIAN_NOISE Adds Gaussian noise to an image.
    %   noisyImage = gaussian_noise(I, snr) adds zero-mean Gaussian noise to 
    %   the input image I with standard deviation calculated as Imax/snr.
    %
    %   noisyImage = gaussian_noise(I, snr, bias) adds Gaussian noise with
    %   the specified mean (bias) and standard deviation Imax/snr.
    %
    %   Inputs:
    %       I    - Input image
    %       snr  - Signal-to-Noise Ratio: SNR = Imax/std  
    %              Higher values = less noise (typical values: 5-100)
    %       bias - Mean value of the Gaussian noise (default: 0)
    %
    %   Output:
    %       noisyImage - Image with added Gaussian noise
    
    % Input validation
    if nargin < 2
        error('Usage: gaussian_noise(image, snr, [bias])');
    end
    
    % Set default bias to 0 if not provided
    if nargin < 3
        bias = 0;
    end
    
    if ~isnumeric(I)
        error('Image must be numeric');
    end
    
    if ~isscalar(snr) || snr <= 0
        error('SNR must be a positive scalar');
    end
    
    if ~isscalar(bias)
        error('Bias must be a scalar value');
    end
    
    % Store original data type
    originalType = class(I);
    
    % Convert to double for processing
    I = double(I);
    
    % % Get image properties
    % maxIntensity = max(I(:));
    % fprintf('Maximum image intensity: %f\n', maxIntensity);
    
    % % Calculate noise standard deviation using SNR = Imax/std formula
    % noiseStd = maxIntensity / snr;
    
    % % Generate Gaussian noise with specified mean and calculated variance
    % noise = bias + noiseStd * randn(size(I));
    % If you want a power-based SNR:
    signal_power = mean(I(:).^2);
    noise_power = signal_power / snr;
    noiseStd = sqrt(noise_power);
    noise = + bias+ noiseStd * randn(size(I));

    % Add noise to the image
    noisyImage = I + noise;
    noisyImage = noisyImage;
    % Limit output values to the range of the original data type
    switch originalType
        case 'uint8'
            noisyImage = max(0, min(255, noisyImage));
          
        case 'uint16'
            noisyImage = max(0, min(65535, noisyImage));
        case 'int16'
            noisyImage = max(-32768, min(32767, noisyImage));
        case 'single'
            % No clipping needed, but ensure it's in the 0-1 range 
            % if the original was in that range
            if max(I(:)) <= 1
                noisyImage = max(0, min(1, noisyImage));
            end
        case 'double'
            % No clipping needed, but ensure it's in the 0-1 range 
            % if the original was in that range
            if max(I(:)) <= 1
                noisyImage = max(0, min(1, noisyImage));
            end
    end
    

end