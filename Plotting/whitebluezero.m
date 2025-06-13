function c = whitebluezero(lower_limit,upper_limit)
% whitebluezero(lower_limit, upper_limit)
% 
% This function generates a colormap that transitions from white to shades of blue.
% The generated colormap is ideal for visualizing data where low values should be 
% represented by white and high values by blue, with a smooth gradient in between.
% It is particularly useful for data that ranges from zero to positive values, where 
% the midpoint (zero) is represented by white, and the positive values are represented 
% in varying intensities of blue.
%
%
% Inputs:
%   lower_limit - The lower limit of the data range, representing the minimum 
%                 value to be visualized in the colormap (typically near zero).
%   upper_limit - The upper limit of the data range, representing the maximum 
%                 value to be visualized in the colormap (typically positive).
%
% Outputs:
%   c           - An M-by-3 matrix representing the RGB values of the colormap.
%                 The colormap transitions from white (for lower values) to blue 
%                 (for upper values) in a smooth gradient.
%
% Description:
%   - The function first calculates the ratio of the lower limit to the upper 
%     limit, which determines how quickly the color transitions from white to blue.
%   - The colormap is created by defining separate color channels for red (r), green 
%     (g), and blue (b). The red and green channels both decrease linearly from the 
%     ratio to zero, while the blue channel remains constant at 1 (blue).
%   - The function returns the resulting colormap matrix, which can be used for 
%     visualizing data with values between the lower and upper limits.





m=256;
ratio=1-(lower_limit/upper_limit);
% From [1 1 1] to [0 0 1]
r = linspace(ratio, 0, m)';
g = linspace(ratio, 0, m)';
b = linspace(1, 1, m)';

c = [r, g, b];
