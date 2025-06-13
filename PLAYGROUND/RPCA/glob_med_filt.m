function [A,B]= glob_med_filt(A,B,sig_bnds)


% glob_med_filt
% Takes the data matrix x and returns a matrix filtered_data of the same
% size with global outliers outside of +-sig_threshold replaced by
% neighbour interpolated data.
% Input
% A            X vectors, array(float)
% B            Y vectors, array(float)
% sig_bnds     number of standard deviations to form a bound, str
% outputs:
% X and Y vector arrays with global outliers removed
 
        
        
    
        % Pre dimension
        x_med = zeros(size(A));
        x_sig = zeros(size(A));
        globfilt = zeros([size(A),2]);
    
        % Calculate the global filter
        
        % Calculate median and standard deviation in the vectorfields
        
        x_med= median(A, 'all');
        x_sig = std(A, 0, 'all'); 
        y_med= median(B, 'all');
        y_sig = std(B, 0, 'all'); 
        
        % Calculate global median filter. If any value is outside of the +- 
        % sig_bnds range in any of the three components, the whole 3-component 
        % vector needs to be replaced.
        globfilt(:, :, 1) =    (A < (x_med - sig_bnds*x_sig) | A > (x_med + sig_bnds*x_sig));
        globfilt(:, :, 2) =    (B < (y_med - sig_bnds*y_sig) | B > (y_med + sig_bnds*y_sig));
        
            

   
    
        globfilt = any(globfilt, 3);
    
   
        % Replace identified filter positions with NaNs.
        x_temp = A;
        x_temp(globfilt) = NaN;
        y_temp = B;
        y_temp(globfilt) = NaN;
        
        % Interpolate infill of NaNs from neighbours. See doc(inpaint_nans) for
        % details.
        A = inpaint_nans(double(x_temp), 3);
        B = inpaint_nans(double(y_temp), 3);
    
       
    
        

end