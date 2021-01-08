% SP_ERROR_VEL_DARCY: Evaluate the error for the velocity for Darcy
%
%   [err, errl2, errl2_bd] = sp_error_vel_darcy (space, msh, u, uex, bnd_sides, solver)
%
% INPUT:
%
%    space:     object defining the space of discrete functions (see sp_vector)
%    msh:       object defining the domain partition and the quadrature rule (see msh_cartesian)
%    u:         vector of dof weights
%    uex:       function handle to evaluate the exact solution
%    bnd_sides: boundaries weakly imposed essential b.c.
%    solver:    "1" or "2" 
%
%
% OUTPUT:
%
%     err:       error in a suitable norm
%
%
% Copyright (C) 2020 Riccardo Puppi
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.

function [err, err_l2, err_l2_bd] = sp_error_vel_darcy (space, msh, u, uex, bnd_sides, solver)


p = space.scalar_spaces{1}.degree(1);

% element size h
msh_elem = msh_evaluate_element_list(msh, 1);
h_el = msh_elem.element_size(1);

 err_l2 = 0;
  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (space, msh_col, 'value', true, 'divergence', true);
    err_l2 = err_l2 +  (sp_l2_error (sp_col, msh_col, u, uex)).^2;
  end
   
  err_l2_bd = 0;
  for iside = bnd_sides
      msh_side = msh_eval_boundary_side(msh, iside);
      msh_side_from_interior = msh_boundary_side_from_interior(msh, iside);
      sp_bnd = space.constructor (msh_side_from_interior);
      sp_bnd = struct (sp_precompute (sp_bnd, msh_side_from_interior, 'value', true));
      if solver==1
        err_l2_bd = err_l2_bd + h_el^(-1)*(sp_error_udotn (sp_bnd, msh_side, u, uex)).^2;
      elseif solver==2
        err_l2_bd = err_l2_bd + h_el^(-p-1)*(sp_error_udotn (sp_bnd, msh_side, u, uex)).^2;
      end
  end
 
err = sqrt(err_l2+err_l2_bd);
err_l2 = sqrt(err_l2);
err_l2_bd = sqrt(err_l2_bd);

end
