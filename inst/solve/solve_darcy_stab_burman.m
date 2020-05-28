% Modified by Riccardo in order to test stabilization of pressures for
% Darcy

% use the norm proposed by Erik Burman for the pressures

% SOLVE_DARCY_STAB_BURMAN: Solve a Darcy problem with a B-spline discretization. 
%
% The function solves the Stokes problem
%
%   k*vel + grad(press) = f    in Omega
%              div(vel) = g    in Omega
%                     p = p_D  on Gamma_N (natural b.c.)
%           vel\cdot \n = U_N  on Gamma_D (essential b.c.)
%
% USAGE:
%
%  [geometry, msh, space_v, vel, space_p, press] = ...
%                        solve_darcy_stab (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with natural boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - f:            force term
%    - g:            mass term
%    - p_D:          function for natural condition (if nmnn_sides is not empty)
%    - u_N:          function for Dirichlet boundary condition
%    - permeability: permeability coefficient (mu in the equation)
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
%
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

function [geometry, msh, space_v, vel, space_p, press, int_dofs, p_dofs, M_ur, M_pr, Ar, B0r, B1r] = ...
                          solve_darcy_stab_burman (problem_data, method_data)

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
[space_v, space_p] = sp_bspline_fluid_mod (element_name, ...
                geometry.nurbs.knots, nsub, degree, regularity, msh);

% Assemble the matrices
if (msh.rdim == 2)
  fun_one = @(x, y) ones (size(x));
elseif (msh.rdim == 3)
  fun_one = @(x, y, z) ones (size(x));
end
A = op_u_v_tp (space_v, space_v, msh, permeability); 
B0 = op_divv_q_tp (space_v, space_p, msh);
E = op_f_v_tp (space_p, msh, fun_one).';
F = op_f_v_tp (space_v, msh, f);
G0 = op_f_v_tp (space_p, msh, g);
% matrix inducing scalar product for pressures
M_p = op_u_v_tp (space_p, space_p, msh);
M_p = op_gradu_gradv_el (space_p, space_p, msh);




% matrix inducing scalar product for velocities
M_u = op_u_v_tp (space_v, space_v, msh) + op_divu_divv_tp (space_v, space_v, msh); 

vel   = zeros (space_v.ndof, 1);
press = zeros (space_p.ndof, 1);

% Apply natural boundary conditions
rhs_nmnn = zeros(space_v.ndof,1);
for iside = nmnn_sides
    msh_side = msh_eval_boundary_side (msh, iside);
    msh_side_from_interior = msh_boundary_side_from_interior (msh, iside);
    sp_side_v = space_v.constructor (msh_side_from_interior);
    sp_side_v = struct (sp_precompute (sp_side_v, msh_side_from_interior, 'value', true));
    sp_side_v.dofs = 1:sp_side_v.ndof;
    x = cell (msh_side.rdim, 1);
    for idim = 1:msh_side.rdim
        x{idim} = squeeze (msh_side.geo_map(idim,:,:));
    end
    pval = reshape (p_D (x{:}, iside), 1, msh_side.nqn, msh_side.nel);
    rhs_nmnn(sp_side_v.dofs) = rhs_nmnn(sp_side_v.dofs) + op_f_vdotn(sp_side_v, msh_side, pval);
end

% Apply essential boundary conditions in a strong sense
if (strcmpi (element_name, 'RT') || strcmpi (element_name, 'NDL'))
  [vel_drchlt, drchlt_dofs] = sp_drchlt_l2_proj_udotn (space_v, msh, drchlt_sides, velex);
  else
  [vel_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space_v, msh, u_N, drchlt_sides);
end

An = spalloc (space_v.ndof, space_v.ndof, 3*space_v.ndof);
Bn = spalloc (space_p.ndof, space_v.ndof, 3*space_v.ndof);
Fn = spalloc (space_v.ndof, 1, space_v.ndof);
Gn = spalloc (space_p.ndof, 1, space_p.ndof);
M_un = spalloc (space_v.ndof, space_v.ndof, 3*space_v.ndof);
M_pn = spalloc (space_p.ndof, space_p.ndof, 3*space_p.ndof);

% Apply essential boundary conditions in a weak sense
for iside = weak_drchlt_sides
    msh_side = msh_eval_boundary_side (msh, iside);
    msh_side_from_interior = msh_boundary_side_from_interior (msh, iside);
    sp_side_v = space_v.constructor (msh_side_from_interior);
    sp_side_v = struct (sp_precompute (sp_side_v, msh_side_from_interior, 'value', true));
    sp_side_p = space_p.constructor (msh_side_from_interior);
    sp_side_p = struct (sp_precompute (sp_side_p, msh_side_from_interior, 'value', true));
    for idim = 1:msh.rdim
        x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
    end
    mat_udotn_vdotn = op_udotn_vdotn (sp_side_v, sp_side_v, msh_side, ones(msh_side.nqn, msh_side.nel));
    An = An + Cpen*mat_udotn_vdotn;
    Bn = Bn + op_q_v_n (sp_side_v, sp_side_p, msh_side);
    vel_times_coeff = bsxfun (@times, velex(x{:}, iside), ones(msh_side.rdim, msh_side.nqn, msh_side.nel));
    Fn = Fn + Cpen*op_fdotn_vdotn (sp_side_v, msh_side, vel_times_coeff);
    Gn = Gn + op_f_v (sp_side_p, msh_side, u_N(x{:}, iside));
    M_un = M_un + mat_udotn_vdotn;
    M_pn = M_pn + op_u_v (sp_side_p, sp_side_p, msh_side, msh_side.charlen);
end
A = A + An;
F = F + Fn;

B1 = B0 - Bn;
G1 = G0 - Gn;
M_u = M_u + M_un;
M_p = M_p + M_pn;


p_dofs = 1:space_p.ndof;
vel(drchlt_dofs) = vel_drchlt;
int_dofs = setdiff (1:space_v.ndof, drchlt_dofs);
nintdofs = numel (int_dofs);
if ~isempty(drchlt_dofs)
    if symmetric
        rhs_dir  = [-A(int_dofs, drchlt_dofs); ...
           B1(p_dofs, drchlt_dofs)]*vel(drchlt_dofs);
    else
        rhs_dir  = [-A(int_dofs, drchlt_dofs); ...
           B0(p_dofs, drchlt_dofs)]*vel(drchlt_dofs);
    end
end

% assemble the final matrix
if symmetric
    mat = [ A(int_dofs, int_dofs), -B1(:,int_dofs).';
        -B1(:,int_dofs),  sparse(size (B1,1), size (B1,1))];
    rhs = [F(int_dofs) - rhs_nmnn(int_dofs);
        -G1 ];
else
    mat = [ A(int_dofs, int_dofs), -B1(:,int_dofs).';
        -B1(:,int_dofs),  sparse(size (B1,1), size (B1,1))];
    rhs = [F(int_dofs) - rhs_nmnn(int_dofs);
        -G0 ];
end

% add strong Dirichlet b.c.
if ~isempty(drchlt_sides)
   rhs = rhs + rhs_dir; 
end

% filter out constant pressures in the case of no natural b.c. 
if isempty (nmnn_sides)
   avg_pres_vec = [spalloc(1, nintdofs, 0), E(p_dofs)];
   mat = [mat, avg_pres_vec'; ...
       avg_pres_vec, 0];
   rhs = [rhs; 0];
end

% solve the linear system
sol = mat \ rhs;
vel(int_dofs) = sol(1:nintdofs);
if isempty (nmnn_sides)
    press = sol(1+nintdofs:end-1);
else
    press = sol(1+nintdofs:end);
end

% restrict to the right dofs
Ar = A(int_dofs, int_dofs);
B0r = B0(p_dofs,int_dofs);
B1r = B1(p_dofs,int_dofs);
M_ur = M_u(int_dofs,int_dofs);
M_pr = M_p(p_dofs,p_dofs);


end

% 
% if method_data.stab
%     % JUST FOR DEBUGGING: matrice di restrizione da non stabilizzato a stabilizzato
%     C = speye(space_p.ndof);
%     % a single bad element
%     % decide bad elements and good neighbours (not on the boundary)
%     bad_el_id = msh.nel/2 - msh.nel_dir(1)/2; 
%     good_nghb_el_id = bad_el_id + msh.nel_dir(1)+1; % north-east (top-right diagonally)
%     C(bad_el_id,bad_el_id)=0;
%     C(good_nghb_el_id,bad_el_id)=1;
%     p_dofs(bad_el_id) = [];
% %     % a whole row of bad elements
% %     good_el_id = [msh.nel/2+1 :msh.nel/2 + msh.nel_dir(1)];
% %     bad_el_id = good_el_id - msh.nel_dir(1);
% %     C(bad_el_id,bad_el_id)=0;
% %     for k = 1:numel(bad_el_id)
% %        C(good_el_id(k), bad_el_id(k)) = 1; 
% %        p_dofs(bad_el_id(k)) = [];
% %     end
%     B = C*B;
%     G = C*G;
%     temp = C*E';
%     E = temp';
%     M_p = C*M_p;
% end

