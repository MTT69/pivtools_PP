function [A,B,median_filter] = loc_med_filt(A,B)

    epsilon_0=0.2;
    eps_threshold=2;
    filt_size=[3 3];
    
    x(:,:,1)=A;
    x(:,:,2)=B;

    [ny, nz] = size(x(:,:,1));
    n_fields = size(x, 3);
    m_filt = filt_size(1); n_filt = filt_size(2);
    % -------------------
    % Local median filter
    % -------------------
    % Calculate local median filter (exclude NaN). Interpolate infill.
    
    for i = 1:n_fields
        % medfilt2 replaces each cell with the meadian of it's m x n
        % neighbours
        x_med(:,:,i)	= medfilt2(x(:,:,i), [m_filt n_filt], 'symmetric');
        
        % Create pad array (array with zeros surrounding the velocity array)
        % Make the pad symmetric and extend in both directions
        pad_size = [floor(m_filt/2) floor(n_filt/2)]; % Expand array with zeros in each direction
        x_pad(:,:,i) = padarray(x(:,:,i), pad_size, 'symmetric', 'both');
    end
    
    % Calculate residual r = u_mn - u_med, excluding center cell
    r_mn = zeros([ny nz m_filt n_filt]);
    
    for m_loop = 1:m_filt
       for n_loop = 1:n_filt
           for i = 1:n_fields
                r_mn(:, :, m_loop, n_loop) = r_mn(:, :, m_loop, n_loop) + (x_med(:,:,i) - x_pad(m_loop-1 + (1:ny), n_loop-1 + (1:nz), i)).^2;
           end
       end
    end
    
    r_mn = reshape(sqrt(r_mn), [ny nz m_filt*n_filt]);
    r_mn = r_mn(:, :, 1:m_filt*n_filt);
    r_med = median(r_mn, 3);
    
    % Norm of the residual between median and unfiltered vector field
    r_2d = zeros(ny, nz);
    for i = 1:n_fields
        r_2d		= r_2d + (x_med(:,:,i) - x(:,:,i)).^2;
    end
    r_2d = sqrt(r_2d);
    
    % Filter then given by threshold
    median_filter	= r_2d ./ (r_med + epsilon_0) > eps_threshold;
    
    % Filter out bad vectors and interpolate to replace
    if any(median_filter(:))

       for i = 1:n_fields
           x_temp = x(:,:, i); 
           x_temp(median_filter) = nan;
           x(:,:,i) = inpaint_nans(double(x_temp));
       end
    
    end

    A=x(:,:,1);
    B=x(:,:,2);

end