% EX_STOKES_LID_DRIVEN_CAVITY_2D_RT: solve the Stokes problem in the unit square (lid-driven cavity) with  Raviart-Thomas elements.
clear all
close all
% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
% problem_data.geo_name = 'geo_square.txt';
problem_data.geo_name = nrbsquare([0,0],1,1,1,0);

% Type of boundary conditions for each side of the domain
problem_data.drchlt_sides = 1:4;
problem_data.nmnn_sides = [];

% Physical parameters
problem_data.viscosity = @(x, y) ones (size (x));

% Force term
fx = @(x, y) zeros(size(x));
fy = @(x, y) zeros(size(x));
problem_data.f  = @(x, y) cat(1, ...
                reshape (fx (x,y), [1, size(x)]), ...
                reshape (fy (x,y), [1, size(x)]));

% Boundary terms
problem_data.h = @(x, y, iside) h(x, y, iside); 


% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
d = 2;
method_data.element_name = 'th';   % Element type for discretization
method_data.degree       = [ d  d];  % Degree of the splines (pressure space)
method_data.regularity   = [ d-1  d-1];  % Regularity of the splines (pressure space)

m = 4;

method_data.nsub         = 2.^[m m];  % Number of subdivisions
method_data.nquad        = [ d+1  d+1];  % Points for the Gaussian quadrature rule

% Penalization parameter for Nitsche's method
factor = 10;
method_data.Cpen = factor*(min(method_data.degree)+1);

% 3) CALL TO THE SOLVER
[geometry, msh, space_v, vel, space_p, press] = ...
                       solve_stokes (problem_data, method_data);

% 4) POST-PROCESSING

% 4.3) PLOT IN MATLAB

% plot velocity solution
figure1 = figure(1);
npts = 20;
vtk_pts = {linspace(0, 1, npts), linspace(0, 1, npts)};
[eu, F] = sp_eval (vel, space_v, geometry, vtk_pts);
[X,  Y] = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
axis equal
title('Computed solution velocity')

% plot pressure solution
figure2 = figure(2);
npts = 128;
vtk_pts = {linspace(0, 1, npts), linspace(0, 1, npts)};
[eu, F] = sp_eval (press, space_p, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
surf (X, Y, eu)
colorbar
view(2)
title('Computed solution pressure')

% 4.2) EXPORT TO PARAVIEW
output_file = 'SQUARE_RT_Deg3_Reg2_Sub10';

fprintf ('The result is saved in the files %s \n and %s \n \n', ...
           [output_file '_vel'], [output_file '_press']);

sp_to_vtk (press, space_p, geometry, vtk_pts, [output_file '_press'], 'press')
sp_to_vtk (vel,   space_v, geometry, vtk_pts, [output_file '_vel'  ], {'velocity', 'divergence'}, {'value', 'divergence'})




function u = h(x, y, iside)
switch(iside)
    case 1
        u = zeros([2, size(x)]);
    case 2
        u = zeros([2, size(x)]);
    case 3
        u1 = -ones(size(x));
        u2 = zeros(size(x));
        u = cat(1, ...
            reshape (u1, [1, size(x)]), ...
            reshape (u2, [1, size(x)]));
    case 4
        u1 = ones(size(x));
        u2 = zeros(size(x));
        u = cat(1, ...
            reshape (u1, [1, size(x)]), ...
            reshape (u2, [1, size(x)]));
end
end