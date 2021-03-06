% OP_F_VDOTN: assemble the vector r = [r(i)], with  r(i) = (f, v_i \cdot n), with n the normal exterior vector.
%
%   rhs = op_f_vdotn (spv, msh, coeff);
%
% INPUT:
%     
%   spv:   structure representing the function space (see sp_vector/sp_eval_boundary_side)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_eval_boundary_side)
%   coeff: source function evaluated at the quadrature points
%
% OUTPUT:
%
%   rhs: assembled right-hand side
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


function rhs = op_f_vdotn (spv, msh, coeff)

shpv = reshape (spv.shape_functions, spv.ncomp, msh.nqn, spv.nsh_max, msh.nel);
coeff = reshape (coeff, 1, msh.nqn, msh.nel);

rhs = zeros (spv.ndof,1);
for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))
        shpv_iel = reshape (shpv(:, :, 1:spv.nsh(iel), iel), spv.ncomp, msh.nqn, spv.nsh(iel), 1);
        shpv_dot_n = sum (bsxfun (@times, shpv_iel, msh.normal(:,:,iel)), 1);
        
        jacdet_weights = reshape (msh.jacdet(:, iel) .* msh.quad_weights(:, iel), 1, msh.nqn);
        coeff_times_jw = bsxfun (@times, jacdet_weights, coeff(:,:,iel));
        
        aux_val = bsxfun (@times, coeff_times_jw, shpv_dot_n);
        %reshape (sum (tmp1, 2), spv.nsh(iel), spu.nsh(iel));
        rhs_loc = sum (sum (aux_val, 1), 2);
        rhs(spv.connectivity(1:spv.nsh(iel), iel)) = rhs(spv.connectivity(1:spv.nsh(iel), iel)) + rhs_loc(:);
    end
end



end

