function [A_hat, E_hat, Y, mu, iter] = inexact_alm_rpc(D, lambda, tol, maxIter,masks,gappy)

% inexact_alm_rpc
%
% This function implements the inexact augmented Lagrange multiplier (ALM) 
% method for Robust Principal Component Analysis (RPCA). It decomposes the 
% input matrix `D` into a low-rank component (`A_hat`) and a sparse error 
% component (`E_hat`) by solving the optimization problem:
%
%   min (A, E) ||A||_* + lambda * ||E||_1 
%       subject to D = A + E
%
% The method can optionally perform selective updates to masked data points
% when the `gappy` mode is enabled, allowing structured or incomplete datasets
% to be handled effectively.
%
% Inputs:
%   D       - m x n matrix of observations/data (required input).
%   lambda  - Regularization parameter controlling the weight of the sparse term
%             in the optimization (default: 1/sqrt(max(m, n))).
%   tol     - Convergence tolerance for the stopping criterion (default: 1e-7).
%   maxIter - Maximum number of iterations allowed (default: 1000).
%   masks   - Logical matrix of size m x n, where `true` indicates uncertain 
%             or missing values to be updated. Used in `gappy` mode.
%   gappy   - Boolean flag indicating if only masked elements should be updated.
%
% Outputs:
%   A_hat   - Low-rank approximation of the input matrix `D`.
%   E_hat   - Sparse error matrix capturing deviations or anomalies in `D`.
%   Y       - Dual variable matrix used in the optimization process.
%   mu      - Final value of the penalty parameter.
%   iter    - Total number of iterations performed until convergence.
%
% Description:
%   - The algorithm iteratively minimizes the augmented Lagrange function:
%
%       L(A, E, Y, mu) = ||A||_* + lambda * ||E||_1 + <Y, D - A - E> + 
%                        (mu / 2) * ||D - A - E||_F^2
%
%   - At each iteration:
%       1. `E_hat` is updated using a soft-thresholding operation.
%       2. `A_hat` is updated using a truncated Singular Value Decomposition 
%          (SVD) on the residuals.
%       3. The dual variable `Y` is updated, and penalty parameter `mu` is scaled.
%
%   - In `gappy` mode, only masked entries in `masks` are updated, leaving 
%     unmasked data unchanged. This is useful for structured reconstruction 
%     or imputation tasks.
%
%   - The process stops when the relative Frobenius norm of the residual 
%     matrix (||D - A_hat - E_hat|| / ||D||) falls below the tolerance `tol`, 
%     or when `maxIter` is reached.
%
% Notes:
%   - The method uses efficient SVD computation (`svdecon`) for better 
%     performance with large matrices.
%   - Progress and diagnostics are printed every 10 iterations, including the 
%     rank of `A_hat`, sparsity of `E_hat`, and convergence metrics.
%
% Example:
%   % Generate sample data
%   D = rand(100, 100) + 5 * eye(100);
%   lambda = 1 / sqrt(size(D, 1));
%   tol = 1e-6;
%   maxIter = 500;
%   masks = rand(100, 100) > 0.8; % Random mask
%   gappy = true;
%
%   % Perform RPCA
%   [A_hat, E_hat, Y, mu, iter] = inexact_alm_rpc(D, lambda, tol, maxIter, masks, gappy);
%
%   % Visualize results
%   imagesc(A_hat); title('Low-Rank Component (A\_hat)');
%
% References:
%   Minming Chen, Microsoft Research Asia (2009). v-minmch@microsoft.com
%   Arvind Ganesh, University of Illinois Urbana-Champaign (abalasu2@illinois.edu)
%   Masking added by Morgan Taylor

%
% Copyright: Perception and Decision Laboratory, University of Illinois, Urbana-Champaign
%            Microsoft Research Asia, Beijing

[m,n] = size(D);

if(nargin < 2) lambda = 1 / sqrt(m); end
if(nargin < 3) tol = 1e-7; elseif(tol == -1) tol = 1e-7; end
if(nargin < 4) maxIter = 1000; elseif(maxIter == -1) maxIter = 1000; end


% initialize
Y = D;
norm_two = norm(Y, 2);
norm_inf = norm( Y(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;

A_hat = zeros( m, n);
E_hat = zeros( m, n);
mu = 1.25/norm_two; % this one can be tuned
mu_bar = mu * 1e7;
rho = 1.5;          % this one can be tuned
d_norm = norm(D, 'fro');

iter = 0;
total_svd = 0;
converged = false;
stopCriterion = 1;
sv = 10;
while ~converged       
    iter = iter + 1;
    
    temp_T = D - A_hat + (1/mu)*Y;
    E_hat = max(temp_T - lambda/mu, 0);
    E_hat = E_hat+min(temp_T + lambda/mu, 0);

%     [U,S,V] = svd(D - E_hat + (1/mu)*Y, 'econ');
    [U,S,V] = svdecon(D - E_hat + (1/mu)*Y,0); % fastest
    
    diagS = diag(S);
    svp = length(find(diagS > 1/mu));
    if svp < sv
        sv = min(svp + 1, n);
    else
        sv = min(svp + round(0.05*n), n);
    end
    
    A_hat = U(:, 1:svp) * diag(diagS(1:svp) - 1/mu) * V(:, 1:svp)';  
    if gappy
        A_hat(~masks) = D(~masks);
    end

    

    total_svd = total_svd + 1;
    
    Z = D - A_hat - E_hat;

%     % Create a figure
%     figure(1)
%     
%     % First subplot
%     subplot(2,2,1); % 2x2 grid, first position
%     imagesc(reshape(A_hat(1:end/2, 1), [149, 319])); 
%     colorbar; % Add colorbar
%     daspect([1 1 1]); % Equal aspect ratio
%     
%     % Second subplot
%     subplot(2,2,2); % 2x2 grid, second position
%     imagesc(reshape(Y(1:end/2, 1), [149, 319])); 
%     colorbar; % Add colorbar
%     daspect([1 1 1]); % Equal aspect ratio
%     
%     % Third subplot
%     subplot(2,2,3); % 2x2 grid, third position
%     imagesc(reshape(E_hat(1:end/2, 1), [149, 319])); 
%     colorbar; % Add colorbar
%     daspect([1 1 1]); % Equal aspect ratio
%     
%     % Fourth subplot
%     subplot(2,2,4); % 2x2 grid, fourth position
%     imagesc(reshape(Z(1:end/2, 1), [149, 319])); 
%     colorbar; % Add colorbar
%     daspect([1 1 1]); % Equal aspect ratio
%     
%     % Adjust the layout (optional)
%     sgtitle('Visualizing Data in Subplots'); % Set a title for the whole figure

    
    
    Y = Y + mu*Z;
    if gappy
        Y(~masks) = 0;
    end
    mu = min(mu*rho, mu_bar);
        
    %% stop Criterion    
    stopCriterion = norm(Z, 'fro') / d_norm;
    if stopCriterion < tol
        converged = true;
    end    
    
    if mod( total_svd, 10) == 0
        disp(['#svd ' num2str(total_svd) ' r(A) ' num2str(rank(A_hat))...
            ' |E|_0 ' num2str(length(find(abs(E_hat)>0)))...
            ' stopCriterion ' num2str(stopCriterion)]);
    end
    
    if ~converged && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
    disp(['#svd ' num2str(total_svd) ' r(A) ' num2str(rank(A_hat))...
            ' |E|_0 ' num2str(length(find(abs(E_hat)>0)))...
            ' stopCriterion ' num2str(stopCriterion)]);
end
