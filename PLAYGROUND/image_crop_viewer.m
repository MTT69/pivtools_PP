
    % IMAGE_CROP_VIEWER Displays the central 128x128 pixels of an input image
    %   IMAGE_CROP_VIEWER(I) takes an image I, extracts the central 128x128
    %   pixels and displays them using imagesc with a range of 0 to 255.
%     I = B00001_A;
    I = filtered_gaussian_frame1;
%     I = filtered_original__frame2;
    % Get the size of the input image
    [height, width, ~] = size(I);
    
    % Check if the image is large enough
    if height < 128 || width < 128
        error('Input image must be at least 128x128 pixels.');
    end
    
    % Calculate the center of the image
    center_h = round(height/2);
    center_w = round(width/2);
    
    % Calculate the indices for the 128x128 central region
    half_size = 64; % Half of 128
    
    row_start = center_h - half_size + 1;
    row_end = center_h + half_size;
    col_start = center_w - half_size + 1;
    col_end = center_w + half_size;
    
    % Extract the central region
    cropped_img = I(row_start:row_end, col_start:col_end, :);
    
    % Display the cropped image with a fixed range [0, 255]
    figure;
    imagesc(cropped_img, [0 255]);
    colormap(gray); % Use grayscale colormap
    colorbar;
    axis image; % Keep the aspect ratio
    title('Central 128x128 Region of Image');
    xlabel('X (pixels)');
    ylabel('Y (pixels)');
