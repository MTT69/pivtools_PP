function p = OMNIPOL(u,v,U,V,dx,mask,nu,rho,num,smooth_dp,smooth_p)
% OMNIPOL - Omnidirectional pressure integration with user-prescribed
% finite number of evenly-distributed polar angles (2D)
%
% This function returns the pressure field using Navier-Stokes equations
% given the instantaneous and mean velocity fields and a specified number
% of polar angles. The unsteady term is approximated using Taylor's
% Hypothesis, therefore the mean velocity field is required as an input.
%
% NOTE: Users must ensure inputs are in appropriate (consistent) units.
%
% Implementation based on 8-path integration code of:
% Dabiri, John O., et al. "An algorithm to estimate unsteady and 
% quasi-steady pressure fields from velocity field measurements." Journal 
% of Experimental Biology 217.3 (2014): 331-336.
% Extension to abritrary number of paths by D. W. Carter (2021)

% Inputs:       u - horizontal velocity field of size [ny nx]
%               v - vertical velocity field of size [ny nx]
%               U - mean horizontal velocity field of size [ny nx]  (use
%                   zeros of size [ny nx] if calculating mean pressure)
%               V - mean vertical velocity field of size [ny nx]  (use
%                   zeros of size [ny nx] if calculating mean pressure)
%              dx - size of (uniform) grid spacing [meters]
%            mask - matrix of ones with NaNs at masked locations [ny nx]
%              nu - kinematic viscosity [meters^2/sec]
%             rho - density [kg / meters^3]
%             num - number of polar angles to consider for integration,
%                   more angles -> slower, less directional bias (default: 8)
%       smooth_dp - option to smooth pressure gradients [true/false]
%        smooth_p - option to smooth pressure [true/false]
%
% Output:       p - instantaneous pressure field

    % Handle defaults
    if nargin < 9
        num = 8;
    end   
    if nargin < 10
        smooth_dp = 1;
    end   
    if nargin < 11
        smooth_p = 1;
    end

    % Obtain matrix dimensions
    [ny,nx] = size(u);

    % First, obtain pressure gradients     
    u_f = u - U;
    v_f = v - V;

    [dudx,dudy] = gradient(u,dx);
    [dvdx,dvdy] = gradient(v,dx);

    [dufdx,dufdy] = gradient(u_f,dx);
    [dvfdx,dvfdy] = gradient(v_f,dx);

    dpdx = -rho*(-(U.*dufdx+V.*dufdy)+(u.*dudx+v.*dudy) ...
        -nu*4.*del2(u,dx,dx)).*mask;
    dpdy = -rho*(-(U.*dvfdx+V.*dvfdy)+(u.*dvdx+v.*dvdy)- ...
        nu*4.*del2(v,dx,dx)).*mask;
    
    % Optional smoothing
    if smooth_dp
        dpdx = smoothNN(dpdx);
        dpdy = smoothNN(dpdy);
    end
    
    % Obtain shifting directions (degrees)
    theta = linspace(0,360,num+1);
    theta = theta(1:end-1);
    
    % Initialize pressure integration storage from each direction
    p_total = zeros(ny,nx,num);
    
    % Integrate pressure gradient along each polar direction
    for i = 1:num
        p_total(:,:,i) = intPolarPGrad(dpdx,dpdy,dx,theta(i));
    end
    
    p = nanmedian(p_total,3).*mask;   
    
    if smooth_p == 1
        p = smoothNN(p).*mask;
    end
    
end

