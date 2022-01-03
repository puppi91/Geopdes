% SOLVE_DARCY: Solve a Darcy problem with a B-spline discretization. 
%
% The function solves the Darcy problem
%
%   k*vel -  grad(press) = f    in Omega
%              div(vel) = g    in Omega
%                     p = p_D  on delta_N (natural b.c.)
%           vel\cdot \n = U_N  on delta_D (essential b.c.)
%
% using the formulation
% 
% (u_h,v_h) + Cpen*<u_h \cdot n, v_h\cdot n>_WE - (p_h, div v_h) 
%  + <p_h, v_h\cdot n>_WE = (f,v_h) + Cpen*<u_N, v_h\cdot n>_WE
%  - <p_D, v_h\cdot n>_N
%  -(q_h, div u_h) + m <q_h, u_h\cdot n> = -(q_h, g) + m <q_h, u_N>
%
% m \in {0,1}: m=0 -> non-symmetric formulation, m=1 -> symmetric
% formulation
%
% USAGE:
%
%  [geometry, msh, space_v, vel, space_p, press] = ...
%                        solve_darcy (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:           name of the file containing the geometry
%    - ntrl_sides:         sides with natural boundary condition (may be empty)
%    - essntl_sides: sides with essential boundary condition
%    - weak_essntl_sides:  sides with essential boundary condition to be
%                          weakly enforced
%    - f:                  force term
%    - g:                  mass term
%    - p_D:                function for natural condition (if ntrl_sides is not empty)
%    - u_N:                function for Dirichlet boundary condition
%    - permeability:       permeability coefficient (mu in the equation)
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

function [geometry, msh, space_v, vel, space_p, press, int_dofs, p_dofs,...
    M_ur, M_pr, Ar, B0r, Br, cond_numb] = ...
                          solve_darcy_prova (problem_data, method_data)

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
  case {'RT_FEM'}
    [~, zeta] = kntrefine (geometry.nurbs.knots, nsub-1, degree, regularity);
end
% Construct msh structure
rule       = msh_gauss_nodes (nquad);
[qn, qw]   = msh_set_quad_nodes (zeta, rule);
msh        = msh_cartesian (zeta, qn, qw, geometry);

% Compute the space structures
[space_v, space_p] = sp_bspline_fluid (element_name, ...
                geometry.nurbs.knots, nsub, degree, regularity, msh);
            
h_el = msh_evaluate_element_list(msh,1).element_size;
p = space_p.degree(1);
% Assemble the matrices
if (msh.rdim == 2)
  fun_one = @(x, y) ones (size(x));
elseif (msh.rdim == 3)
  fun_one = @(x, y, z) ones (size(x));
end
A = op_u_v_tp (space_v, space_v, msh, permeability);
A = A + delta*op_divu_divv_tp (space_v, space_v, msh);
M_u = A;
B = op_divv_p_tp (space_v, space_p, msh);
E = op_f_v_tp (space_p, msh, fun_one).';
F = op_f_v_tp (space_v, msh, f);
G = op_f_v_tp (space_p, msh, divvelex);
G0 = G;
B0 = B;

% matrix inducing scalar product for pressures
M_p = op_gradu_gradv_tp (space_p, space_p, msh);
M_p_jump = zeros(size(M_p));
B_jump = zeros(size(B));
% adding jumps' contribution
for ii = 1:msh.ndim
    for jj = 1: numel(msh.breaks{ii})-2
        [~, msh_left, msh_right] =  msh_on_internal_face (msh, ii, jj+1);
        sp_left_p = space_p.constructor (msh_left);
        sp_left_p = struct (sp_precompute (sp_left_p, msh_left));
        sp_right_p = space_p.constructor (msh_right);
        sp_right_p = struct (sp_precompute (sp_right_p, msh_right));
        sp_left_v = space_v.constructor (msh_left);
        sp_left_v = struct (sp_precompute (sp_left_v, msh_left));
       % msh_face = msh_precompute(msh_left);  % choose one of the two
        msh_face = msh_evaluate_element_list(msh_left, 1:msh_left.nel, 'normal', true, 'idir',ii);
        coeff = ones(msh_face.nqn, msh_face.nel);
        sp_side_p.dofs = 1:sp_left_p.ndof; % choose one of the two
        sp_side_v.dofs = 1:sp_left_v.ndof; % choose one of the two
        M_p_jump(sp_side_p.dofs,sp_side_p.dofs)  = M_p_jump(sp_side_p.dofs,sp_side_p.dofs)  + ...
            + (h_el)^(-1)*(op_u_v (sp_left_p, sp_left_p, msh_face,coeff)...
            - op_u_v (sp_left_p, sp_right_p, msh_face,coeff)...
        - op_u_v (sp_right_p,sp_left_p, msh_face,coeff) + op_u_v (sp_right_p, sp_right_p, msh_face,coeff));
        B_jump(sp_side_p.dofs,sp_side_v.dofs)  = B_jump(sp_side_p.dofs,sp_side_v.dofs)  + ...
            - op_p_vdotn (sp_left_v, sp_left_p, msh_face)...
            + op_p_vdotn (sp_left_v, sp_right_p, msh_face);
    end
end
M_p = M_p + M_p_jump;


vel   = zeros (space_v.ndof, 1);
press = zeros (space_p.ndof, 1);

% Apply natural boundary conditions
rhs_ntrl = zeros(space_v.ndof,1);
for iside = ntrl_sides
    msh_side = msh_eval_boundary_side (msh, iside);
    if (strcmpi (element_name, 'RT_FEM'))
        msh_side_from_interior = msh_boundary_side_from_interior (msh, iside);
        sp_side_v = space_v.constructor (msh_side_from_interior);
        sp_side_v = struct (sp_precompute (sp_side_v, msh_side_from_interior, 'value', true));
        sp_side_v.dofs = 1:sp_side_v.ndof;
    end
    x = cell (msh_side.rdim, 1);
    for idim = 1:msh_side.rdim
        x{idim} = squeeze (msh_side.geo_map(idim,:,:));
    end
    pval = reshape (p_D (x{:}, iside), 1, msh_side.nqn, msh_side.nel);
    rhs_ntrl(sp_side_v.dofs) = rhs_ntrl(sp_side_v.dofs) + op_f_vdotn(sp_side_v, msh_side, pval);
end
F = F + rhs_ntrl;

An = spalloc (space_v.ndof, space_v.ndof, 3*space_v.ndof);
Bn = spalloc (space_p.ndof, space_v.ndof, 3*space_v.ndof);
Fn = spalloc (space_v.ndof, 1, space_v.ndof);
Gn = spalloc (space_p.ndof, 1, space_p.ndof);
M_un = spalloc (space_v.ndof, space_v.ndof, 3*space_v.ndof);
M_pn = spalloc (space_p.ndof, space_p.ndof, 3*space_p.ndof);

% Apply essential boundary conditions in a weak sense
for iside = weak_essntl_sides
    msh_side = msh_eval_boundary_side (msh, iside);
    msh_side_from_interior = msh_boundary_side_from_interior (msh, iside);
    sp_side_v = space_v.constructor (msh_side_from_interior);
    sp_side_v = struct (sp_precompute (sp_side_v, msh_side_from_interior, 'value', true, 'divergence', true));
    sp_side_p = space_p.constructor (msh_side_from_interior);
    sp_side_p = struct (sp_precompute (sp_side_p, msh_side_from_interior, 'value', true));
    for idim = 1:msh.rdim
        x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
    end
    mat_udotn_vdotn = op_udotn_vdotn (sp_side_v, sp_side_v, msh_side, ones(msh_side.nqn, msh_side.nel));
    mat_divu_vdotn = op_divu_vdotn (sp_side_v, sp_side_v, msh_side, ones(msh_side.nqn, msh_side.nel));
    vel_times_coeff = bsxfun (@times, velex(x{:}, iside), ones(msh_side.rdim, msh_side.nqn, msh_side.nel));
    if solver==1
        An = An + Cpen*(h_el^(-1))*mat_udotn_vdotn;
        An = An - delta*mat_divu_vdotn;
        Bn = Bn + op_p_vdotn (sp_side_v, sp_side_p, msh_side);
        Fn = Fn + Cpen*(h_el^(-1))*op_fdotn_vdotn (sp_side_v, msh_side, vel_times_coeff);
        Gn = Gn + op_fdotn_v (sp_side_p, msh_side, velex(x{:}, iside));
        M_un = M_un + Cpen*(h_el^(-1))*mat_udotn_vdotn;
    elseif solver==2
        An = An + Cpen*(h_el^(-p-1))*mat_udotn_vdotn;
        An = An - delta*mat_divu_vdotn;
        Fn = Fn + Cpen*(h_el^(-p-1))*op_fdotn_vdotn (sp_side_v, msh_side, vel_times_coeff);%...
             +(~isempty(ntrl_sides))*op_f_vdotn(sp_side_v, msh_side, pressex(x{:}));
        M_un = M_un + Cpen*(h_el^(-p-1))*mat_udotn_vdotn;
        M_pn = M_pn + (h_el^(-1))*op_u_v (sp_side_p, sp_side_p, msh_side, msh_side.charlen);
    end
end
A = A + An;
F = F + Fn;
B = B - Bn;
if (~exist('symmetric','var'))
    symmetric = true;
end
if symmetric
   G = G - Gn; 
end
M_u = M_u + M_un;
M_p = M_p + M_pn;


%------------------------------DEBUGGING----------------------------------
% non_trimmed_elem_ids = 1:msh.nel;
% trimmed_elems = {};
% B_jump_bis = op_vdotn_qjump_trimming_bis (space_v,space_p, msh, non_trimmed_elem_ids, trimmed_elems);
% max(max(abs(B_jump - B_jump_bis)))
% 
% M_p_jump_bis = h_el^(-1)*op_ujump_vjump_trimming_bis (space_p, msh, non_trimmed_elem_ids, trimmed_elems);
% max(max(abs(M_p_jump - M_p_jump_bis)))

fprintf('----------------check i.b.p----------------');
B_grad = op_v_gradp_tp (space_v, space_p, msh);
B_tilde = -B_grad + B_jump;
max(max(abs(B-B_tilde)))

% check that jumps are zero when pressures are continuous
fprintf('----------------check pressures jumps----------------');
max(max(abs(M_p_jump)))
max(max(abs(B_jump)))

%--------------------------------------------------------------------------


% Apply essential boundary conditions in a strong sense
if (strcmpi (element_name, 'RT_FEM'))
  [vel_essntl, essntl_dofs] = sp_drchlt_l2_proj_udotn (space_v, msh, essntl_sides, velex);
  else
  [vel_essntl, essntl_dofs] = sp_drchlt_l2_proj (space_v, msh, u_N, essntl_sides);
end

p_dofs = 1:space_p.ndof;
vel(essntl_dofs) = vel_essntl;
int_dofs = setdiff (1:space_v.ndof, essntl_dofs);
nintdofs = numel (int_dofs);

if ~isempty(essntl_dofs)
    if symmetric
        rhs_dir  = [-A(int_dofs, essntl_dofs); ...
           -B(:, essntl_dofs)]*vel(essntl_dofs);
    else
        rhs_dir  = [-A(int_dofs, essntl_dofs); ...
           -B0(:, essntl_dofs)]*vel(essntl_dofs);
    end
end

% assemble linear system
Ar = A(int_dofs, int_dofs);
B0r = B0(:,int_dofs);
Br = B(:,int_dofs);
M_ur = M_u(int_dofs,int_dofs);
M_pr = M_p;
if~isempty(weak_essntl_sides)
   Bnr = Bn(:, int_dofs);
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
if ~symmetric
    mat = [ Ar, Br.';
        B0r, spalloc(space_p.ndof, space_p.ndof, 0)];
    rhs = [Fr;
        G0r];
else
    mat = [ Ar, Br.';
        Br, spalloc(space_p.ndof, space_p.ndof, 0)];
    rhs = [Fr;
        Gr];
end

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
if ~isempty(ntrl_sides)
    press = sol(1+nintdofs : end);
else
    press= sol(1+nintdofs :end-1);
end

% condition number global matrix
cond_numb = condest(mat);
% err = sp_error_press_darcy (space_p, msh, press, pressex, gradpressex, weak_essntl_sides)

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

