% OP_G_DIVV: assemble the right-hand side vector r = [r(i)], with  r(i) = (g, div v_i).
%
%   rhs = op_g_divv (spv, msh, coeff);
%
% INPUT:
%
%   spv:   structure representing the function space (see sp_vector/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   coeff: source function evaluated at the quadrature points
%
% OUTPUT:
%
%   rhs: assembled right-hand side
% 
% Copyright (C) 2021 Riccardo Puppi
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

function rhs = op_g_divv (spv, msh, coeff)
  
 coeff = reshape (coeff, 1, msh.nqn, msh.nel);

 rhs   = zeros (spv.ndof, 1);
 shpv  = reshape (spv.shape_function_divs, 1, msh.nqn, spv.nsh_max, msh.nel);

 for iel = 1:msh.nel
   if (all (msh.jacdet(:,iel)))
     jacdet_weights = reshape (msh.jacdet(:, iel) .* msh.quad_weights(:, iel), 1, msh.nqn);

     coeff_times_jw = bsxfun (@times, jacdet_weights, coeff(:,:,iel));

     shpv_iel = reshape (shpv(:, :, 1:spv.nsh(iel), iel), 1, msh.nqn, spv.nsh(iel));

     aux_val = bsxfun (@times, coeff_times_jw, shpv_iel);
     rhs_loc = sum (sum (aux_val, 1), 2);
     rhs(spv.connectivity(1:spv.nsh(iel), iel)) = rhs(spv.connectivity(1:spv.nsh(iel), iel)) + rhs_loc(:); 
   else
     warning ('geopdes:jacdet_zero_at_quad_node', 'op_g_divv: singular map in element number %d', iel)
   end
 end
 
end


