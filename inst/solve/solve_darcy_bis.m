% SOLVE_DARCY_BIS: Solve a Darcy problem with a B-spline discretization. 
%
% The function solves the Stokes problem
%
%   k*vel -  grad(press) = f    in Omega
%              div(vel) = g    in Omega
%                     p = p_D  on Gamma_N (natural b.c.)
%           vel\cdot \n = U_N  on Gamma_D (essential b.c.)
%
% using the formulation
% 
% (u_h,v_h) + Cpen*<u_h \cdot n, v_h\cdot n>_WE - (p_h, div v_h) 
%   = (f,v_h) + <Cpen*u_N+p_D, v_h\cdot n>_WE 
%  - <p_D, v_h\cdot n>_N
%  -(q_h, div u_h)  = -(q_h, g) 
%
%
% USAGE:
%
%  [geometry, msh, space_v, vel, space_p, press] = ...
%                        solve_darcy_bis (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - ntrl_sides:   sides with natural boundary condition (may be empty)
%    - essntl_sides: sides with essential boundary condition
%    - f:            force term
%    - g:            mass term
%    - p_D:          function for natural condition (if ntrl_sides is not empty)
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
%     - Cpen:        penalty parameter (for the weak imposition of the essential b.c.)
%
% OUTPUT:
%
%  geometry:  geometry structure (see geo_load)
%  msh:       mesh object that defines the quadrature rule (see msh_cartesian)
%  space_v:   space object for the velocity (see sp_vector)
%  vel:       the computed degrees of freedom for the velocity
%  space_p:   space object for the pressure (see sp_scalar)
%  press:     the computed degrees of freedom for the pressure
%  cond_numb: condition number
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

function [geometry, msh, space_v, vel, space_p, press, int_dofs, p_dofs, M_ur, M_pr, Ar, B0r, cond_numb] = ...
                          solve_darcy_bis (problem_data, method_data)

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
  case {'RT', 'RT_FEM', 'TH', 'NDL'}
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
A = op_u_v_tp (space_v, space_v, msh, permeability); 
B0 = op_divv_q_tp (space_v, space_p, msh);
E = op_f_v_tp (space_p, msh, fun_one).';
F = op_f_v_tp (space_v, msh, f);
G0 = op_f_v_tp (space_p, msh, divvelex);

h_el = msh_evaluate_element_list(msh,1).element_size;
% matrix inducing scalar product for pressures
M_p = op_gradu_gradv_tp (space_p, space_p, msh);

for ii = 1:msh.ndim
    for jj = 1: numel(msh.breaks{ii})-2
        [~, msh_left, msh_right] =  msh_on_internal_face (msh, ii, jj+1);
        sp_left = space_p.constructor (msh_left);
        sp_left = struct (sp_precompute (sp_left, msh_left));
        sp_right = space_p.constructor (msh_right);
        sp_right = struct (sp_precompute (sp_right, msh_right));
        msh_face = msh_precompute(msh_left);  % choose one of the two
        coeff = ones(msh_face.nqn, msh_face.nel);
        sp_side_p.dofs = 1:sp_left.ndof; % choose one of the two
        M_p(sp_side_p.dofs,sp_side_p.dofs) = M_p(sp_side_p.dofs,sp_side_p.dofs)...
            + (h_el)^(-1)*(op_u_v (sp_left, sp_left, msh_face,coeff) - op_u_v (sp_left, sp_right, msh_face,coeff)...
        - op_u_v (sp_right,sp_left, msh_face,coeff) + op_u_v (sp_right, sp_right, msh_face,coeff));     
    end
end

for iside = weak_essntl_sides
    msh_side = msh_eval_boundary_side (msh, iside);
    sp_side_p = space_p.constructor (msh_side);
    sp_side_p = struct (sp_precompute (sp_side_p, msh_side));
    sp_side_p.dofs = 1:sp_side_p.ndof;
    coeff = ones(msh_side.nqn, msh_side.nel);
    M_p(sp_side_p.dofs,sp_side_p.dofs) =  M_p(sp_side_p.dofs,sp_side_p.dofs)...
        +  (h_el)^(-1-method_data.degree(1))*op_u_v (sp_side_p, sp_side_p, msh_side,coeff);
end

% matrix inducing scalar product for velocities
M_u = op_u_v_tp (space_v, space_v, msh); 
for iside = weak_essntl_sides
    msh_side = msh_eval_boundary_side (msh, iside);
    msh_side_from_interior = msh_boundary_side_from_interior (msh, iside);
    sp_side_v = space_v.constructor (msh_side_from_interior);
    sp_side_v = struct (sp_precompute (sp_side_v, msh_side_from_interior, 'value', true));
    sp_side_v.dofs = 1:sp_side_v.ndof;
    coeff = ones(msh_side.nqn, msh_side.nel);
    M_u(sp_side_v.dofs,sp_side_v.dofs) =  M_u(sp_side_v.dofs,sp_side_v.dofs)...
        +  (h_el)^(-1-method_data.degree(1))*op_udotn_vdotn (sp_side_v, sp_side_v, msh_side,coeff); 
end

vel   = zeros (space_v.ndof, 1);
press = zeros (space_p.ndof, 1);

% Apply natural boundary conditions
rhs_ntrl = zeros(space_v.ndof,1);
for iside = ntrl_sides
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
    rhs_ntrl(sp_side_v.dofs) = rhs_ntrl(sp_side_v.dofs) + op_f_vdotn(sp_side_v, msh_side, pval);
end

% Apply essential boundary conditions in a strong sense
if (strcmpi (element_name, 'RT') || (strcmpi (element_name, 'RT_FEM')) || strcmpi (element_name, 'NDL'))
  [vel_essntl, essntl_dofs] = sp_drchlt_l2_proj_udotn (space_v, msh, essntl_sides, velex);
  else
  [vel_essntl, essntl_dofs] = sp_drchlt_l2_proj (space_v, msh, u_N, essntl_sides);
end

An = spalloc (space_v.ndof, space_v.ndof, 3*space_v.ndof);
Fn = spalloc (space_v.ndof, 1, space_v.ndof);
M_un = spalloc (space_v.ndof, space_v.ndof, 3*space_v.ndof);
M_pn = spalloc (space_p.ndof, space_p.ndof, 3*space_p.ndof);

% Apply essential boundary conditions in a weak sense
for iside = weak_essntl_sides
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
    vel_times_coeff = bsxfun (@times, velex(x{:}, iside), ones(msh_side.rdim, msh_side.nqn, msh_side.nel));
    pval = reshape (p_D (x{:}, iside), 1, msh_side.nqn, msh_side.nel);
    Fn = Fn + Cpen*op_fdotn_vdotn (sp_side_v, msh_side, vel_times_coeff)...
        + op_f_vdotn(sp_side_v, msh_side, pval);
    M_un = M_un + mat_udotn_vdotn;
    M_pn = M_pn + op_u_v (sp_side_p, sp_side_p, msh_side, msh_side.charlen);
end
A = A + An;
F = F + Fn;

M_u = M_u + M_un;
M_p = M_p + M_pn;

p_dofs = 1:space_p.ndof;
vel(essntl_dofs) = vel_essntl;
int_dofs = setdiff (1:space_v.ndof, essntl_dofs);
nintdofs = numel (int_dofs);
if ~isempty(essntl_dofs)
    rhs_dir  = [-A(int_dofs, essntl_dofs); ...
        -B0(p_dofs, essntl_dofs)]*vel(essntl_dofs);
end

% assemble the final matrix

% restrict to the right dofs
Ar = A(int_dofs, int_dofs);
B0r = B0(p_dofs,int_dofs);
M_ur = M_u(int_dofs,int_dofs);
M_pr = M_p(p_dofs,p_dofs);
Er = E(p_dofs);

Ar = 0.5*(Ar+Ar'); % symmetrize
mat = [ Ar, B0r.';
    B0r,  sparse(space_p.ndof, space_p.ndof)];
rhs = [F(int_dofs) + rhs_ntrl(int_dofs);
    G0 ];

% add strong essential b.c.
if ~isempty(essntl_sides)
   rhs = rhs + rhs_dir; 
end

% filter out constant pressures in the case of no natural b.c. 
if isempty (ntrl_sides)
   avg_pres_vec = [spalloc(1, nintdofs, 0), Er];
   mat = [mat, avg_pres_vec'; ...
       avg_pres_vec, 0];
   rhs = [rhs; 0];
end

% solve the linear system
sol = mat \ rhs;
vel(int_dofs) = sol(1:nintdofs);
if isempty (ntrl_sides)
    press = sol(1+nintdofs:end-1);
else
    press = sol(1+nintdofs:end);
end

% condition number global matrix
cond_numb = condest(mat);

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

