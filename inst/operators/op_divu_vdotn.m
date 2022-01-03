% OP_DIVU_VDOTN: assemble the matrix M = [m(i,j)], m(i,j) = (div u_j, v_i * n), with n the exterior normal vector.
%
%   mat = op_divu_vdotn (spu, spv, msh, coeff);
%   [rows, cols, values] = op_divu_vdotn (spu, spv, msh, coeff);
%
% INPUT:
%
%  spu:   structure representing the space of trial functions (see sp_vector/sp_eval_boundary_side)
%  spv:   structure representing the space of test functions  (see sp_vector/sp_eval_boundary_side)
%  msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_eval_boundary_side)
%  coeff: physical parameter
%
% OUTPUT:
%
%   mat:    assembled matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
%
% Copyright (C) 2021    Riccardo Puppi
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

function varargout = op_divu_vdotn (spu, spv, msh, coeff)

divshpu = reshape (spu.shape_function_divs, msh.nqn, spu.nsh_max, msh.nel);
shpv = reshape (spv.shape_functions, spv.ncomp, msh.nqn, spv.nsh_max, msh.nel);

rows = zeros (msh.nel * spv.nsh_max * spu.nsh_max, 1);
cols = zeros (msh.nel * spv.nsh_max * spu.nsh_max, 1);
values = zeros (msh.nel * spv.nsh_max * spu.nsh_max, 1);

jacdet_weights = msh.jacdet .* msh.quad_weights .* coeff;

ncounter = 0;
for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
        divshpu_iel = reshape (divshpu(:, 1:spu.nsh(iel), iel), msh.nqn, spu.nsh(iel), 1);
        shpv_iel = reshape (shpv(:, :, 1:spv.nsh(iel), iel), spv.ncomp, msh.nqn, spv.nsh(iel), 1);
        n_iel = reshape(msh.normal(:,:, iel), msh.rdim, msh.nqn, 1);
        
        v_cdot_n = reshape(sum(bsxfun(@times, shpv_iel, n_iel), 1), msh.nqn, spv.nsh(iel), 1);
        
        jacdet_iel = reshape(jacdet_weights(:,iel), msh.nqn, 1, 1);
        
        v_cdot_n_times_j =  reshape(bsxfun(@times, v_cdot_n, jacdet_iel), msh.nqn, 1, spv.nsh(iel));
        
        tmp = sum(bsxfun(@times, v_cdot_n_times_j, divshpu_iel), 1);
        
        elementary_values = reshape (tmp, spu.nsh(iel), spv.nsh(iel));
        
        [rows_loc, cols_loc] = ndgrid (spu.connectivity(:,iel), spv.connectivity(:,iel));
        indices = rows_loc & cols_loc;
        rows(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = rows_loc(indices);
        cols(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = cols_loc(indices);
        values(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = elementary_values(indices);
        ncounter = ncounter + spu.nsh(iel)*spv.nsh(iel);
        
    else
       warning ('geopdes:jacdet_zero_at_quad_node', 'op_divu_vdotn: singular map in element number %d', iel)
    end
end

if (nargout == 1 || nargout == 0)
    varargout{1} = sparse (rows, cols, values, spu.ndof, spv.ndof);
elseif (nargout == 3)
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = values;
else
    error ('op_divu_vdotn: wrong number of output arguments')
end


end

