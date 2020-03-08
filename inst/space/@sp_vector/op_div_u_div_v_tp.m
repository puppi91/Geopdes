% OP_DIV_U_DIV_V_TP: assemble the matrix B = [b(i,j)], b(i,j) = (div v_j, div v_i), exploiting the tensor product structure.
%
%   mat = op_div_u_div_v_tp (spu, spv, msh);
%   [rows, cols, values] = op_div_u_div_v_tp (spu, spv, msh);
%
% INPUT:
%
%   spu:     object representing the space of trial functions (see sp_vector)
%   spv:     object representing the space of test functions (see sp_vector)
%   msh:     object defining the domain partition and the quadrature rule (see msh_cartesian)
%
% OUTPUT:
%
%   mat:    assembled matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
% Copyright (C) 2020, Riccardo Puppi
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

function varargout = op_div_u_div_v_tp (space1, space2, msh)

A = spalloc (space2.ndof, space1.ndof, 3*space1.ndof);

for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp1_col = sp_evaluate_col (space1, msh_col,'value', false, 'divergence', true);
    sp2_col = sp_evaluate_col (space2, msh_col, 'value', false, 'divergence', true);
    
    A = A + op_div_u_div_v (sp1_col, sp2_col, msh_col);
end

if (nargout == 1)
    varargout{1} = A;
elseif (nargout == 3)
    [rows, cols, vals] = find (A);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
end

end
