set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(groot,'defaultTextInterpreter','latex')
clear variables
close all
clc
addpath(genpath('aux_matlab'));
addpath(genpath('BaseFlows'));
verbose         = true;
mean_fields = load(fullfile('C:\Users\mtt1e23\OneDrive - University of Southampton\Documents\#current_processing\PINNs\90\Processed\ENSEMBLE\Resolvent','resolvent_mean.mat'));

% Physical parameters
baseFlow.Re     = 3.8e4;          % Reynolds number
baseFlow.Ma     = 0.03;          % Mach number
baseFlow.Pr     = 0.7;          % Prandtl number
baseFlow.kappa 	= 1.4;          % heat capacity ratio
baseFlow.T_0 	= 273.15 ;      % temperature


nEig            = 3;            % Arnoldi method number of eigenvalues

%% Coordinate axes
nx              = 300;
ny              = 100;
dx              = 0.05;
x_min           = 0;
x_max           = 15;
y_min           = 0;
y_i             = 0.05;     
y_max           = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create mesh and obtain differentiation matrices                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTES: -why have CreateMesh and DeformMesh separately?
%        -put y_symmetry,x_periodicity,alpha_filter in options and default
%         to false?
%        -missing default for sponge function
%        -opts for eigs hardwired in resolvent function

% Cartesian mesh in computational domain
y_symmetry      = false;        % use symmetry on y coordinate around y=0 (for axysymmetric problems)
x_periodicity   = false;        % use periodic b.c. on x
alpha_filter    = 0;            % spatial filter coefficient
xrange          = [x_min x_max];% domain range in x
yrange          = [y_min y_max];       % domain range in y
FDorder         = 4;
mesh           = CreateMesh(xrange,yrange,nx,ny,FDorder, ...     
                             y_symmetry,x_periodicity,alpha_filter); %construct mesh                     

x               = mesh.X;           
y               = mesh.Y;
mask = (y >= 0 & y <= 1) & (x >= 0 & x <= 3);
mesh            = MeshMask(mesh,mask);

%% Blasius Solution

baseFlow.U = mean_fields.ux_large(1:5:end, 1:5:end);
baseFlow.U(1,:)=0;
baseFlow.V = mean_fields.uy_large(1:5:end, 1:5:end);
baseFlow.V(1,:)=0;
nan_mask = isnan(baseFlow.U);
baseFlow.RHO = ones(size(baseFlow.U));
baseFlow.RHO(nan_mask)=NaN;
baseFlow.W = zeros(size(baseFlow.U));
baseFlow.W(nan_mask)=NaN;
baseFlow.T = ones(size(baseFlow.U));
baseFlow.T(nan_mask)=NaN;

baseFlow        = sutherland_air(baseFlow);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sponge                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sponge defined in the origial, non-deformed, mesh.
x = mesh.X;
y = mesh.Y;

xs_trans    = .2;  xs_offset =  .2;
ys_trans    = 1 ; ys_offset =  1;
spongeAmp   = 5;
ind = mesh.usedInd;

mesh.sponge = nan(size(mesh.X));
mesh.sponge(ind) = ...
        spongeAmp/2*max(tanh( ( y(ind)-max(y(:)) + ys_offset )/ys_trans)+1, ...
                        tanh(-( x(ind)-min(x(:)) - xs_offset)/xs_trans)+1 + ...
                        tanh( ( x(ind)-max(x(:)) + xs_offset)/xs_trans)+1  );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualize base flow                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
    figure('name','Base Flow')
    vars = {baseFlow.RHO ,'$\rho$';
            baseFlow.U,'$U$';
            baseFlow.V,'$V$';
            baseFlow.W,'$W$';
            baseFlow.T,'$T$';           
            baseFlow.MU,'$\mu$';
            mesh.sponge,'Sponge';};
    plotFlow(mesh.X,mesh.Y,vars,4,2);
    drawnow
end
if verbose
    figure
   
    title('Physical domain')
    hold on
    ind = mesh.usedInd;
    plot(mesh.X(ind), mesh.Y(ind), '.k','HandleVisibility','off');
    
    p=mesh.usedInd;
    for bondaries= fields(mesh.idx)'
        ids = mesh.idx.(bondaries{1});
        plot(mesh.X(p(ids)), mesh.Y(p(ids)),'o');
    end
    legend(fields(mesh.idx),'Location', 'Best')
    xlabel('$x$');
    ylabel('$y$');
    axis equal tight
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up linear operator and boundary conditions                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    dqdt = L0 q
%%
[L0,idx] 	= GetLinProblem(mesh,baseFlow,'2D',0);

% Enforce Dirichlet b.c.s for u,v,w and T on all boundaries
borders='lrbt';  vars = 'uvwT';
[L0,idx_dirchlet] = BC_Dirichlet(L0,idx,borders,vars);

[W,invW] = GetCompEnergyNorm(mesh,baseFlow,'2D');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Resolvent analysis                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dq/dt = L*q + B*v
% u = C*q
% v: input/forcing, u:output/response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input/forcing matrix
B                       = ones(mesh.ngp*5,1);
B([idx.T_j;idx.rho_j;idx.w_j])  = 0; % no temperature and density forcing
B( idx_dirchlet  )              = 0; % no forcing on Dirichlet b.c. 
B                               = spdiags(B,0,mesh.ngp*5,mesh.ngp*5);

% Output/observable matrix
C                       = B;    

% Compute optimal forcings and responses with and without filters (for
% comparison)

% Define omega values in two segments
% omega_values_initial = 1:10; % Omega values from 1 to 10
% omega_values_later = 10:10:150; % Omega values from 10 to 150 in steps of 10
% 
% % Concatenate the two arrays
% omega_values = [omega_values_initial, omega_values_later];

omega_values = 1.1:0.1:2.9;


% Define the variable BFS (assuming it's already defined in your script)
BFS_ = 'BFS_90'; % Define BFS according to your specific use case

% Initialize array to store maximum gains
max_gains = zeros(size(omega_values)); % Preallocate for efficiency

% Define other required variables (L0, nEig, W, invW, B, C, mesh.filters, mesh.X, mesh.Y, etc.)
% These should be set up before running the loop

for i = 1:length(omega_values)
    omega = omega_values(i);
    
    % Perform resolvent calculation
    [gains, responses, forces] = resolvent(L0, omega, nEig, W, invW, B, C, mesh.filters);

    % Store the maximum gain for this omega
    max_gains(i) = max(gains);
    
    % Generate file name based on BFS and omega
    file_name = sprintf('%s_omega_%d.mat', BFS_, omega);
    
    % Save the variables to a .mat file
    save(file_name, 'gains', 'responses', 'forces');
    
    % Plot mode gains
    if verbose
        % Create full screen figure for mode gains
        fig1 = figure('Name', 'Mode gains', 'NumberTitle', 'off');
        set(fig1, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); % Full screen
        
        bar(gains);
        xlabel('Mode');
        ylabel('Gain');
        title('Resolvent gains');
        
        % Save plots
        saveas(fig1, sprintf('%s_mode_gains_omega_%d.fig', BFS_, omega));
        saveas(fig1, sprintf('%s_mode_gains_omega_%d.jpg', BFS_, omega));
        
        % Plot resolvent forcing and response modes for u_j
        fig2 = figure('Name', 'Resolvent forcing and response modes for u_j', 'NumberTitle', 'off');
        set(fig2, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); % Full screen
        
        vars_u = {real(forces(idx.u_j,1)), '$f_u^{(1)}$'; real(responses(idx.u_j,1)), '$u^{(1)}$';
                  real(forces(idx.u_j,2)), '$f_u^{(2)}$'; real(responses(idx.u_j,2)), '$u^{(2)}$';
                  real(forces(idx.u_j,3)), '$f_u^{(3)}$'; real(responses(idx.u_j,3)), '$u^{(3)}$'};
        
        axs_u = plotFlow(mesh.X, mesh.Y, vars_u, 3, 2, mesh.usedInd);
        
        % Set axis limits for u_j plots
        for j = 1:length(axs_u)
            ylim(axs_u(j), [y_min, y_max]); % Set y scale for all plots
            xlim(axs_u(j), [x_min, x_max]); % Set x scale for all plots
        end
        
        % Save plots
        saveas(fig2, sprintf('%s_u_j_modes_omega_%d.fig', BFS_, omega));
        saveas(fig2, sprintf('%s_u_j_modes_omega_%d.jpg', BFS_, omega));
        
        % Plot resolvent forcing and response modes for v_j
        fig3 = figure('Name', 'Resolvent forcing and response modes for v_j', 'NumberTitle', 'off');
        set(fig3, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); % Full screen
        
        vars_v = {real(forces(idx.v_j,1)), '$f_v^{(1)}$'; real(responses(idx.v_j,1)), '$v^{(1)}$';
                  real(forces(idx.v_j,2)), '$f_v^{(2)}$'; real(responses(idx.v_j,2)), '$v^{(2)}$';
                  real(forces(idx.v_j,3)), '$f_v^{(3)}$'; real(responses(idx.v_j,3)), '$v^{(3)}$'};
        
        axs_v = plotFlow(mesh.X, mesh.Y, vars_v, 3, 2, mesh.usedInd);
        
        % Set axis limits for v_j plots
        for j = 1:length(axs_v)
            ylim(axs_v(j), [y_min, y_max]); % Set y scale for all plots
            xlim(axs_v(j), [x_min, x_max]); % Set x scale for all plots
        end
        
        % Save plots
        saveas(fig3, sprintf('%s_v_j_modes_omega_%d.fig', BFS_, omega));
        saveas(fig3, sprintf('%s_v_j_modes_omega_%d.jpg', BFS_, omega));
        
        % Close figures
        close(fig1);
        close(fig2);
        close(fig3);
    end
end

% Plot and save the maximum gains against omega
fig4 = figure('Name', 'Maximum Gain vs Omega', 'NumberTitle', 'off');
set(fig4, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); % Full screen

plot(omega_values, max_gains, '-o');
xlabel('Omega');
ylabel('Maximum Gain');
title('Maximum Gain vs Omega');

% Save the plot
saveas(fig4, sprintf('%s_max_gains_vs_omega.fig', BFS_));
saveas(fig4, sprintf('%s_max_gains_vs_omega.jpg', BFS_));

% Close the figure
close(fig4);



