function lambda_ci = swirling(u, v, dx, dy)
    % swirling_strength calculates the swirling strength (lambda_ci)
    % for a 2D velocity field (u,v) given grid spacings dx and dy.
    %
    % NOTE: dx and dy should be computed from your coordinate system.
    % For example:
    %   dx = gradient(Co_ords.Co_ords(i).x);
    %   [~, dy] = gradient(Co_ords.Co_ords(i).y);
    %   dy = -dy;  % Correct for your coordinate system
    %
    % Syntax:
    %   lambda_ci = swirling_strength(u, v, dx, dy)
    %
    % Inputs:
    %   u, v  - 2D matrices of the velocity components.
    %   dx    - Grid spacing in the x-direction.
    %   dy    - Grid spacing in the y-direction (after correction).
    %
    % Output:
    %   lambda_ci - 2D matrix of swirling strength values.
    
    % Set default grid spacing if not provided
    if nargin < 3, dx = 1; end
    if nargin < 4, dy = 1; end
    
    % Compute spatial gradients using the provided grid spacings
    [du_dx, du_dy] = gradient(u, dx(1), dy(1));
    [dv_dx, dv_dy] = gradient(v, dx(1), dy(1));

    du_dy = -du_dy;
    dv_dy = -dv_dy;
    
    % Compute the trace and determinant:
    trace_val = du_dx + dv_dy;
    det_val   = du_dx .* dv_dy - du_dy .* dv_dx;
    
    % Compute the discriminant of the characteristic equation
    discriminant = trace_val.^2 - 4 * det_val;
    
    % For regions with swirling motion, the discriminant is negative.
    % The swirling strength is defined as:
    %   lambda_ci = 0.5 * sqrt( -discriminant ) for discriminant < 0.
    % Here, we ensure non-negative arguments using max(0,-discriminant).
    lambda_ci = 0.5 * sqrt( max(0, -discriminant) );
end
