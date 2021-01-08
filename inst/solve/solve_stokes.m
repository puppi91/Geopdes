% SOLVE_STOKES: Solve a Stokes problem with a B-spline discretization. 
%
% The function solves the Stokes problem
%
%   -div(mu(x) grad(vel)) + grad(press) = f    in Omega
%                              div(vel) = 0    in Omega
%             mu(x) dvel/dn - press * n = g    on Gamma_N
%                                   vel = h    on Gamma_D
%
% USAGE:
%
%  [geometry, msh, space_v, vel, space_p, press] = ...
%                        solve_stokes (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with natural boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - f:            force term
%    - g:            function for natural condition (if nmnn_sides is not empty)
%    - h:            function for Dirichlet boundary condition
%    - viscosity:    viscosity coefficient (mu in the equation)
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:       degree of the spline functions for pressure
%    - regularity:   continuity of the spline functions for pressure
%    - nsub:       number of subelements with respect to the geometry mesh 
%                   for the pressure space (nsub=1 leaves the mesh unchanged)
%    - nquad:        number of points for Gaussian quadrature rule
%    - element_name: one of {TH,SG,RT,NDL}, specify how to build the velocity
%                    space from the data for the pressure space
%                     +TH  is the generalized Taylor-Hood element
%                     +SG  is the SubGrid element
%                     +RT  the generalized Raviart-Thomas element
%                     +NDL the generalized Nedelec element of the second class
%                    for more details see the references below
% OUTPUT:
%
%  geometry: geometry structure (see geo_load)
%  msh:      mesh object that defines the quadrature rule (see msh_cartesian)
%  space_v:  space object for the velocity (see sp_vector)
%  vel:      the computed degrees of freedom for the velocity
%  space_p:  space object for the pressure (see sp_scalar)
%  press:    the computed degrees of freedom for the pressure
%
% See also EX_STOKES_SQUARE_* EX_STOKES_ANNULUS_* for some examples.
%
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
% Copyright (C) 2011, Andrea Bressan, Rafael Vazquez
% Copyright (C) 2014, 2015  Rafael Vazquez
% Copyright (C) 2020 Riccardo Puppi
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [geometry, msh, space_v, vel, space_p, press, drchlt_dofs] = ...
                          solve_stokes (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% load geometry
geometry = geo_load (geo_name);

% Compute the mesh structure using the finest mesh
switch (upper(element_name))
  case {'RT', 'TH', 'NDL'}
    [~, zeta] = kntrefine (geometry.nurbs.knots, nsub-1, degree, regularity);
  case {'SG'}
    [~, zeta] = kntrefine (geometry.nurbs.knots, 2*nsub-1, degree, regularity);
end
rule       = msh_gauss_nodes (nquad);
[qn, qw]   = msh_set_quad_nodes (zeta, rule);
msh        = msh_cartesian (zeta, qn, qw, geometry);

% Compute the space structures
[space_v, space_p] = sp_bspline_fluid (element_name, ...
                geometry.nurbs.knots, nsub, degree, regularity, msh);

% Assemble the matrices
if (msh.rdim == 2)
  fun_one = @(x, y) ones (size(x));
elseif (msh.rdim == 3)
  fun_one = @(x, y, z) ones (size(x));
end
A = op_gradu_gradv_tp (space_v, space_v, msh, viscosity); 
M_u = A;
B = op_divv_q_tp (space_v, space_p, msh);
M_p = op_u_v_tp (space_p, space_p, msh);
G = op_f_v_tp (space_p, msh, divvelex);
E = op_f_v_tp (space_p, msh, fun_one).';
G0 = G;
B0 = B;
F = op_f_v_tp (space_v, msh, f);

vel   = zeros (space_v.ndof, 1);
press = zeros (space_p.ndof, 1);

% Apply natural boundary conditions
rhs_nmnn = zeros(space_v.ndof,1);
for iside = nmnn_sides
  msh_side = msh_eval_boundary_side (msh, iside);
  if (strcmpi (element_name, 'RT') || strcmpi (element_name, 'NDL'))
    msh_side_from_interior = msh_boundary_side_from_interior (msh, iside);

    sp_side = space_v.constructor (msh_side_from_interior);
    sp_side = struct (sp_precompute (sp_side, msh_side_from_interior, 'value', true));
    sp_side.dofs = 1:sp_side.ndof;
  else
    sp_side  = sp_eval_boundary_side (space_v, msh_side);
  end
  x = cell (msh_side.rdim, 1);
  for idim = 1:msh_side.rdim
    x{idim} = squeeze (msh_side.geo_map(idim,:,:));
  end
  gval = reshape (g (x{:}, iside), msh.rdim, msh_side.nqn, msh_side.nel);

  rhs_nmnn(sp_side.dofs) = rhs_nmnn(sp_side.dofs) + op_f_v (sp_side, msh_side, gval);
end
F = F + rhs_nmnn;

% Apply essential boundary conditions in a weak sense
An = spalloc (space_v.ndof, space_v.ndof, 3*space_v.ndof);
Bn = spalloc (space_p.ndof, space_v.ndof, 3*space_v.ndof);
Fn = spalloc (space_v.ndof, 1, space_v.ndof);
Gn = spalloc (space_p.ndof, 1, space_p.ndof);
M_un = spalloc (space_v.ndof, space_v.ndof, 3*space_v.ndof);
M_pn = spalloc (space_p.ndof, space_p.ndof, 3*space_p.ndof);
for iside = weak_drchlt_sides
    msh_side = msh_eval_boundary_side (msh, iside);
    msh_side_from_interior = msh_boundary_side_from_interior (msh, iside);
    sp_side_v = space_v.constructor (msh_side_from_interior);
    sp_side_v = struct (sp_precompute (sp_side_v, msh_side_from_interior, 'value', true,'gradient',true));
    sp_side_p = space_p.constructor (msh_side_from_interior);
    sp_side_p = struct (sp_precompute (sp_side_p, msh_side_from_interior, 'value', true));
    for idim = 1:msh.rdim
        x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
    end
    coeff_at_qn = viscosity (x{:});
    % assemble A_n, M_un, F_n
    Bbd = op_gradv_n_u (sp_side_v, sp_side_v, msh_side, coeff_at_qn);
    g_times_coeff = bsxfun (@times, h(x{:}, iside), ...
        reshape(coeff_at_qn ,[1, msh_side.nqn, msh_side.nel]));
    gradv_n_g = op_gradv_n_f (sp_side_v, msh_side, g_times_coeff);
    coeff_at_qn =  coeff_at_qn * Cpen ./ msh_side.charlen;    
    C = op_u_v (sp_side_v, sp_side_v, msh_side, coeff_at_qn);
    g_times_coeff = bsxfun (@times, h(x{:}, iside), ...
        reshape(coeff_at_qn ,[1, msh_side.nqn, msh_side.nel]));
    g_cdot_v = op_f_v (sp_side_v, msh_side, g_times_coeff);
    An = An + (Bbd + Bbd' - C);
    M_un = M_un + op_u_v (sp_side_v, sp_side_v, msh_side, 1./msh_side.charlen);
    Fn = Fn - gradv_n_g + g_cdot_v;
    
    % assemble B_n, M_qn, G_n
    Bn = Bn + op_q_v_n(sp_side_v, sp_side_p, msh_side);
    Gn = Gn + op_fdotn_v(sp_side_p, msh_side, h (x{:}, iside));
    M_pn = M_pn + op_u_v(sp_side_p, sp_side_p, msh_side, msh_side.charlen);
end
A = A - An;
F = F + Fn;
B = B - Bn;
if (~exist ('symmetric', 'var'))
    symmetric = true;
end
if symmetric
    G = G - Gn;
end
M_u = M_u + M_un;
M_p = M_p + M_pn; 


% Apply essential  boundary conditions in a strong sense. For RT and NDL elements the normal
%  component is imposed strongly, and the tangential one is imposed weakly.
if (strcmpi (element_name, 'RT') || strcmpi (element_name, 'NDL'))
  [A_mat_wD, A_rhs_wD, M_u_add] = ...
      my_sp_weak_drchlt_bc_stokes (space_v, msh, drchlt_sides, h, viscosity, Cpen);
  A = A - A_mat_wD;
  F = F + A_rhs_wD;
  M_u = M_u + M_u_add;
  [vel_drchlt, drchlt_dofs] = sp_drchlt_l2_proj_udotn (space_v, msh, drchlt_sides, h);
else
  [vel_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space_v, msh, h, drchlt_sides);
end

vel(drchlt_dofs) = vel_drchlt;
int_dofs = setdiff (1:space_v.ndof, drchlt_dofs);
nintdofs = numel (int_dofs);

if symmetric
    rhs_dir  = [-A(int_dofs, drchlt_dofs);...
        B(:, drchlt_dofs)] * vel(drchlt_dofs);
else
    rhs_dir  = [-A(int_dofs, drchlt_dofs); ...
        B0(:, drchlt_dofs)] * vel(drchlt_dofs);
end

% assemble linear system
Ar = A(int_dofs, int_dofs);
Br = B(:, int_dofs);
B0r = B0(:, int_dofs); % volumetric term

M_ur = M_u(int_dofs, int_dofs);
M_pr = M_p;

if ~isempty(weak_drchlt_sides)
    Bnr = Bn(:, int_dofs); % boundary term
    Gnr = Gn;
else
    Bnr = zeros(size(Br));
    Gnr = zeros(space_p.ndof,1);
end
Gr = G;
G0r = G0;
Fr = F(int_dofs);
Er = E;

Ar = 0.5*(Ar+Ar'); % symmetrize
if(~symmetric)
    mat = [ Ar, -Br.';
        -B0r, spalloc(space_p.ndof, space_p.ndof, 0)];
    rhs = [Fr;
        -G0r];
else
    mat = [ Ar, -Br.';
        -Br, spalloc(space_p.ndof, space_p.ndof, 0)];
    rhs = [Fr;
        -Gr];
end

% add strong essential b.c.
if ~isempty(drchlt_sides)
    rhs = rhs + rhs_dir;
end

% filter out constant pressures in the case of no natural b.c. 
if isempty(nmnn_sides) 
    avg_pres_vec = [spalloc(1, nintdofs, 0), Er];
    mat = [mat, avg_pres_vec'; ...
        avg_pres_vec, 0];
    rhs = [rhs; 0];
end

% solve the linear system
sol = mat \ rhs;
vel(int_dofs) = sol(1:nintdofs);
if ~isempty(nmnn_sides)
    press = sol(1+nintdofs : end);
else
    press= sol(1+nintdofs :end-1);
end

end
