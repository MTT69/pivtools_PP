function P = EU_FDM(X, Y, U, V, U_MEAN, V_MEAN, dUdt, dVdt, rho, nu)
% This script does pressure reconstruction from a planar particle image
% velocimetry (PIV) field. It first obtains the time-averaged pressure
% gradients from statistics of PIV measurements by means of the
% Reynolds-averaged Navier-Stokes equation, as detailed in van Oudheusden
% et al. (2007).
% EU_FDM: 

% INPUT: Time-averaged velocity fields and Reynold's stresses in the
% streamwise and wall-normal directions. Coordinate vectors must be sorted
% according to the vector fields.
% OUTPUT: It returns a .MAT file of the reconstructed pressure field.

% M. FERREIRA
% January 17, 2018
% Version 2.0 (succeeding version, v1.3)
% - Each solver has a dedicated function;
% - Gradient function grad2D: forward/backward differences of the boundary 
% nodes + central differences on the interior nodes;
% - Function nodeID includes MESH function from previous versions;
% D W CARTER
% October 5, 2021
% Added unsteady terms instead of invoking Taylor's Hypothesis.

%% USER DEFINED VARIABLES --------------------------------------------------
dlet_BC = 0;        % 0. N-BC   1. Apply D-BC
% dlet_ind   = 	{sub2ind(size(U_MEAN)+1, n:m, p:q), [], [], ...};	% specify (n:m) & (p,q)
% dlet_d     = 	{'north', [], [], ...};                             % boundary
% dlet_value =	{0*ones(1,261), [], [], ...};                       % pressure values
% if further boundaries are to be selected, keep adding to the cell array

%% MAIN PROGRAMME ----------------------------------------------------------
% step 1: build nodes and connectivity matrices
pad = 1;
U_MEAN = padarray(U_MEAN,[pad,pad],nan,'both');
V_MEAN = padarray(V_MEAN,[pad,pad],nan,'both');
dUdt = padarray(dUdt,[pad,pad],nan,'both');
dVdt = padarray(dVdt,[pad,pad],nan,'both');
U = padarray(U,[pad,pad],nan,'both');
V = padarray(V,[pad,pad],nan,'both');
h = abs(X(2)-X(1));
% + Identify type of boundary nodes
node = nodeID(U_MEAN);
[dP_dx, dP_dy, f] = EU(X,Y,U_MEAN,V_MEAN,U,V,dUdt,dVdt,rho,nu,node);

% step 2: Apply boundary conditions
% + Unless specified, Neumann conditions are applied at the
% boundary by default. Dirichlet BC are stated below **USER INPUT**
[a2,a1] = size(U_MEAN);
if dlet_BC == 1
    dlet_ind = {sub2ind(size(f),(a2-1)*ones(1,a1-4),3:a1-2)};	% nodes
    dlet_d = {'south'};                                 % direction
    dlet_value = {0.5*rho*( mean(U_MEAN(end-1,3:a1-2))^2  - (U(end-1,3:a1-2).^2 + V(end-1,3:a1-2).^2) )};
    [dlet_c, dlet_n] = dletBC(dlet_ind, dlet_d, dlet_value, node);
else
    dlet_c = [];
    dlet_n = [];
end

% step 3: Builds load vector and stiffness matrix
% + load vector
b = build_b(node, dlet_c, dlet_n, dP_dx, dP_dy, f, h);
id = 1:length(b);

% + stiffness matrix
[I,J,V] = build_A(node, dlet_c, dlet_n);
clearvars -except X Y U_MEAN pad dlet_BC f b I J V id
A = sparse(I,J,V, max(I), max(J),numel(V));

% + the original indexing of the nodes changed due to the padding
% of the vector fields, so the forcing vector and stiffness matrix
% contain blank rows and columns corresponding to the nan values
% removed above.
aux3 = find(isnan(f(:))| f(:)==0);
aux3 = aux3(aux3<length(b));
A(aux3, :) = [];  %rows
A(:, aux3) = [];  %columns
b(aux3) = [];
id(aux3) = [];

% step 4: Solve system of equations
p = A\b';

% + reconstruct pressure field
[idr,idc] = ind2sub(size(f),id);
P = full(sparse(idr,idc,p));
P(P==0) = nan;
if dlet_BC == 0; P = P-nanmean(nanmean(P)); end

% + unpad/pad arrays (revert back to original size)
P = padarray(P,size(U_MEAN)-size(P),nan,'post');
P = P(pad+1:end-pad,pad+1:end-pad);
