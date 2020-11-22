% EX_STOKES_SQUARE_RT: solve the Stokes problem in the unit square with generalized Raviart-Thomas elements.
% modified by Riccardo in order to have convergence error rates
clear all
close all
% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
% problem_data.geo_name = 'geo_square.txt';
problem_data.geo_name = nrb4surf([0,0], [1,0], [0,1], [1,2]);

% Type of boundary conditions for each side of the domain
problem_data.drchlt_sides = [];
problem_data.nmnn_sides = [1];
problem_data.weak_drchlt_sides = [2:4];

% Physical parameters
problem_data.viscosity = @(x, y) ones (size (x));

% Force term
% fx = @(x, y) (6 * x + y .* cos(x .* y) + 2 * cos(y) .* sin(x));
% fy = @(x, y) (x .* cos(x .* y) - 2 * cos(x) .* sin(y));
fx = @(x, y) (1 + 2 * cos(y) .* sin(x));
fy = @(x, y) (- 2 * cos(x) .* sin(y));
problem_data.f  = @(x, y) cat(1, ...
                reshape (fx (x,y), [1, size(x)]), ...
                reshape (fy (x,y), [1, size(x)]));
% Boundary terms
uxex = @(x, y) (sin(x) .* cos(y));
uyex = @(x, y) (-sin(y) .* cos(x));

problem_data.h = @(x, y, iside) cat(1, ...
                  reshape (uxex (x,y), [1, size(x)]), ...
                  reshape (uyex (x,y), [1, size(x)]));
% Exact solution, to compute the errors
problem_data.velex = @(x, y) cat(1, ...
                  reshape (uxex (x,y), [1, size(x)]), ...
                  reshape (uyex (x,y), [1, size(x)]));
%problem_data.pressex = @(x, y) (3 * x.^2 + sin(x .* y)) - 1.239811742000564725943866);
problem_data.pressex = @(x, y) x - 5/6*2/3;
problem_data.gradvelex = @test_stokes_square_bc_graduex;
problem_data.divvelex = @(x,y) zeros(size(x));
problem_data.g = @test_stokes_square_neumann_bc; % Neumann bc

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.element_name = 'rt';   % Element type for discretization
method_data.degree       = [ 1  1];  % Degree of the splines (pressure space)
method_data.regularity   = [ 0  0];  % Regularity of the splines (pressure space)
metod_data.symmetric = false;

mrange = [1:6];
i = 1;
for m = mrange

method_data.nsub         = 2.^[m m];  % Number of subdivisions
method_data.nquad        = method_data.degree+1;  % Points for the Gaussian quadrature rule

% Penalization parameter for Nitsche's method
factor = 10;
method_data.Cpen = factor*(min(method_data.degree)+1);

% 3) CALL TO THE SOLVER
[geometry, msh, space_v, vel, space_p, press] = ...
                       solve_stokes (problem_data, method_data);

% 4) POST-PROCESSING
% 4.1) COMPARISON WITH EXACT SOLUTION
error_l2_p(i) = sp_l2_error (space_p, msh, press, problem_data.pressex);
[error_h1_v(i), error_l2_v(i)] = ...
   sp_h1_error (space_v, msh, vel, problem_data.velex, problem_data.gradvelex);
[~,~,error_hdivs_v(i)] = sp_hdiv_error(space_v, msh, vel, problem_data.velex, ...
    problem_data.divvelex);

i = i + 1;
end

% 4.2) EXPORT TO PARAVIEW
output_file = 'SQUARE_RT_Deg3_Reg2_Sub10';

fprintf ('The result is saved in the files %s \n and %s \n \n', ...
           [output_file '_vel'], [output_file '_press']);
vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
sp_to_vtk (press, space_p, geometry, vtk_pts, [output_file '_press'], 'press')
sp_to_vtk (vel,   space_v, geometry, vtk_pts, [output_file '_vel'  ], {'velocity', 'divergence'}, {'value', 'divergence'})

% 4.3) PLOT IN MATLAB
[eu, F] = sp_eval (vel, space_v, geometry, vtk_pts);
[X,  Y] = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));

% plot velocity solution
figure1 = figure(1);
subplot (1,2,2)
eu2 = problem_data.velex (X, Y);
quiver (X, Y, squeeze(eu2(1,:,:)), squeeze(eu2(2,:,:)))
axis equal
title('Exact solution')
subplot (1,2,1)
quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
axis equal
title('Computed velocity')

% plot pressure solution
figure2 = figure(2);
[eu, F] = sp_eval (press, space_p, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
subplot (1,2,2)
eu2 = problem_data.pressex (X, Y);
surf (X, Y, eu2)
axis equal
title('Exact solution')
subplot (1,2,1)
surf (X, Y, eu)
axis equal
title('Computed pressure')

% plot divergence velocity solution
figure3 = figure(3);
[eu, F] = sp_eval (vel, space_v, geometry, vtk_pts,'divergence');
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
subplot (1,2,2)
eu2 = problem_data.divvelex (X, Y);
surf (X, Y, eu2)
axis equal
title('Exact solution')
subplot (1,2,1)
surf (X, Y, eu)
axis equal
title('Computed divergence velocity')


% plot error velocity
figure4 = figure(4);
h = 2.^-mrange;
loglog(h, error_h1_v, 'ko-', h, error_l2_v, 'bo-', h, h.^(method_data.degree(1)),'k--',...
    h, h.^(method_data.degree(1)+1),'b--')
legend('H1 error', 'L2 error', ...
    strcat('h^', num2str(method_data.degree(1))),...
    strcat('h^', num2str(method_data.degree(1)+1)),...
    'Location', 'Best')
orderl2_vel = diff(log(error_l2_v))./diff(log(h));
orderh1_vel = diff(log(error_h1_v))./diff(log(h));
xlabel('h')
title('error plot velocity')
% plot error pressure
figure5 = figure(5);
h = 2.^-mrange;
loglog(h, error_l2_p, 'ko-', ...
    h, h.^(method_data.degree(1)+1),'b--')
legend('L2 error', ...
    strcat('h^', num2str(method_data.degree(1)+1)),...
    'Location', 'Best')
orderl2_press = diff(log(error_l2_p))./diff(log(h));
xlabel('h')
title('error plot pressure')

%!demo
%! ex_stokes_square_rt

%!test
%! problem_data.geo_name = 'geo_square.txt';
%! problem_data.drchlt_sides = 1:4;
%! problem_data.nmnn_sides = [];
%! problem_data.viscosity = @(x, y) ones (size (x));
%! fx = @(x, y) (6 * x + y .* cos(x .* y) + 2 * cos(y) .* sin(x));
%! fy = @(x, y) (x .* cos(x .* y) - 2 * cos(x) .* sin(y));
%! problem_data.f  = @(x, y) cat(1, ...
%!                 reshape (fx (x,y), [1, size(x)]), ...
%!                 reshape (fy (x,y), [1, size(x)]));
%! uxex = @(x, y) (sin(x) .* cos(y));
%! uyex = @(x, y) (-sin(y) .* cos(x));
%! problem_data.h = @(x, y, iside) cat(1, ...
%!                   reshape (uxex (x,y), [1, size(x)]), ...
%!                   reshape (uyex (x,y), [1, size(x)]));
%! problem_data.velex = @(x, y) cat(1, ...
%!                   reshape (uxex (x,y), [1, size(x)]), ...
%!                   reshape (uyex (x,y), [1, size(x)]));
%! problem_data.pressex = @(x, y) (3 * x.^2 + sin(x .* y) - 1.239811742000564725943866);
%! problem_data.gradvelex = @test_stokes_square_bc_graduex;
%! method_data.element_name = 'rt';   % Element type for discretization
%! method_data.degree       = [ 3  3];  % Degree of the splines (pressure space)
%! method_data.regularity   = [ 2  2];  % Regularity of the splines (pressure space)
%! method_data.nsub         = [10 10];  % Number of subdivisions
%! method_data.nquad        = [ 5  5];  % Points for the Gaussian quadrature rule
%! factor = 10;
%! method_data.Cpen = factor*(min(method_data.degree)+1);
%! [geometry, msh, space_v, vel, space_p, press] = ...
%!                        solve_stokes (problem_data, method_data);
%! error_l2_p = sp_l2_error (space_p, msh, press, problem_data.pressex);
%! [error_h1_v, error_l2_v] = ...
%!    sp_h1_error (space_v, msh, vel, problem_data.velex, problem_data.gradvelex);
%! div = sp_eval (vel, space_v, geometry, [20 20], 'divergence');
%! assert (msh.nel, 100)
%! assert (space_p.ndof, 169)
%! assert (space_v.ndof, 364)
%! assert (error_l2_p, 4.93866840636771e-08, 1e-15)
%! assert (error_h1_v, 3.50559200757476e-06, 1e-15)
%! assert (error_l2_v, 5.53697801153355e-08, 1e-15)
%! assert (max (abs (div(:))) < 1e-12)