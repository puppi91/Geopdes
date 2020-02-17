% EX_LAPLACE_SQUARE: solve the Poisson problem in the unit square with a B-spline discretization.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_square.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4];

% Physical parameters
problem_data.c_diff  = @(x, y) ones(size(x));

x0 = 0.5;
y0 = 0.5;
a = 0.4;
b = 0.4;

% Source and boundary terms
problem_data.f = @(x, y) f(x, y);
problem_data.h = @(x, y, ind) h(x, y, ind);

% Exact solution (optional)
problem_data.uex     = @(x, y) uex(x, y);
problem_data.graduex = @(x, y) graduex(x,y);


d = 4;
mrange = [1:10];

i = 1;
for m = mrange
    % 2) CHOICE OF THE DISCRETIZATION PARAMETERS
    clear method_data
    method_data.degree     = [d d];       % Degree of the splines
    method_data.regularity = [d-1 d-1];       % Regularity of the splines
    method_data.nsub       = [m m];       % Number of subdivisions
    method_data.nquad      = [d+1 d+1];       % Points for the Gaussian quadrature rule
    
    % 3) CALL TO THE SOLVER
    
    [geometry, msh, space, u] = solve_laplace (problem_data, method_data);
    
    
    % Display errors of the computed solution in the L2 and H1 norm
    [error_h1(i), error_l2(i)] = ...
        sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex);
    i = i + 1;
end


% 4.2) PLOT IN MATLAB. COMPARISON WITH THE EXACT SOLUTION

% 4) POST-PROCESSING
% 4.1) EXPORT TO PARAVIEW

output_file = 'Square_BSP_Deg3_Reg2_Sub9';
vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')
figure1 = figure(1);
[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
subplot (1,2,1)
surf (X, Y, eu)
title ('Numerical solution'), axis tight
subplot (1,2,2)
surf (X, Y, problem_data.uex (X,Y))
title ('Exact solution'), axis tight


% plot error velocity
figure2 = figure(2);
hrange = 2.^-mrange;
loglog(hrange, error_h1, 'ko-', hrange, error_l2, 'bo-', hrange, hrange.^(method_data.degree(1)),'k--',...
    hrange, hrange.^(method_data.degree(1)+1),'b--')
legend('H1 error', 'L2 error', ...
    strcat('h^', num2str(method_data.degree(1))),...
    strcat('h^', num2str(method_data.degree(1)+1)),...
    'Location', 'Best')
disp('order L2')
diff(log(error_l2))./diff(log(hrange))
disp('order H1')
diff(log(error_h1))./diff(log(hrange))
xlabel('h')
title('error plot')

function u = uex(x, y)
x0 = 0.5;
y0 = 0.5;
a = 0.4;
b = 0.4;
% t = (x-x0).*(x-x0)./(a.*a)+(y-y0).*(y-y0)./(b.*b);
% u = exp(4./(t-1));
% u = u.*(abs((x-x0)./a).^2+abs((y-y0)./b).^2<1);
u = exp(4./( ((x-x0)/a).^2  + ((y-y0)/b).^2 -1)) .* (abs((x-x0)./a).^2+abs((y-y0)./b).^2<1);
end

function u = h(x, y, ind)
u = zeros(size(x));
end

function u = f(x, y)
x0 = 0.5;
y0 = 0.5;
a = 0.4;
b = 0.4;
% t = (x-x0).*(x-x0)./(a.*a)+(y-y0).*(y-y0)./(b.*b);
% u = exp(4./(t-1));
% x2 = (x-x0).*(x-x0)./(a.*a);
% y2 = (y-y0).*(y-y0)./(b.*b);
% c = 8./(t-1).^4.*u;
% dudx2 = c.*(3.*x2.*x2-(y2-1).*(y2-1)+2.*x2.*(y2+3));
% dudy2 = c.*(3.*y2.*y2-(x2-1).*(x2-1)+2.*y2.*(x2+3));
% uxx = dudx2./(a.*a);
% uyy = dudy2./(b.*b);
% laplacian = - (uxx+uyy);
% u = laplacian.*(abs((x-x0)./a).^2+abs((y-y0)./b).^2<1);
u = (8.*a.^2.*b.^2.*exp(4./( ((x-x0)/a).^2  + ((y-y0)/b).^2 -1)).*(a.^6.*(b.^4-6.*b.^2.*(y-y0).^2-3.*(y-y0).^4)... 
    +a.^4.*(b.^6-2.*b.^4.*(x.^2-2.*x.*x0+x0.^2+(y-y0).^2)...   
        +b.^2.*(y-y0).^2.*(-2.*x.^2+4.*x.*x0-2.*x0.^2+(y-y0).^2))...     
    +a.^2.*b.^4.*(x-x0).^2.*(-6.*b.^2+x.^2-2.*x.*x0+x0.^2-2.*y.^2+4.*y.*y0-2.*y0.^2)...   
    -3.*b.^6.*(x-x0).^4))./(a.^2.*((y-y0).^2-b.^2)+b.^2.*(x-x0).^2).^4;
u = u.* (abs((x-x0)./a).^2+abs((y-y0)./b).^2<1);

end

function u = graduex(x, y)
x0 = 0.5;
y0 = 0.5;
a = 0.4;
b = 0.4;
u =  cat (1, ...
    reshape (-8*(x-x0).*exp(4./( ((x-x0)/a).^2  + ((y-y0)/b).^2 -1))./(a^2*((x-x0).^2./a^2+ (y-y0).^2./b^2 -1).^2).* (abs((x-x0)./a).^2+abs((y-y0)./b).^2<1), [1, size(x)]), ...
    reshape (-8*(y-y0).*exp(4./( ((x-x0)/a).^2  + ((y-y0)/b).^2 -1))./(b^2*((x-x0).^2./a^2+ (y-y0).^2./b^2 -1).^2).* (abs((x-x0)./a).^2+abs((y-y0)./b).^2<1), [1, size(x)]));
end

