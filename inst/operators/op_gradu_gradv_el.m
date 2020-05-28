% OP_GRADU_GRADV_EL: assemble the stiffness matrix A = [a(i,j)], a(i,j) = sum_K (epsilon grad u_j, grad v_i), 
%
%   mat = op_gradu_gradv_el (spu, spv, msh, [epsilon]);
%   [rows, cols, values] = op_gradu_gradv_el (spu, spv, msh, [epsilon]);
%
% INPUT:
%
%   spu:     object representing the space of trial functions (see sp_scalar)
%   spv:     object representing the space of test functions (see sp_scalar)
%   msh:     object defining the domain partition and the quadrature rule (see msh_cartesian)
%   epsilon: function handle to compute the diffusion coefficient (optional)
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix
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

function varargout = op_gradu_gradv_el (space1, space2, msh, coeff)

A = spalloc(space2.ndof, space1.ndof, 3*space1.ndof);

for iel = msh.nel
    msh_el = msh_evaluate_element_list (msh, iel);
    sp1_el = sp_evaluate_element_list (space1, msh_el, 'value', false, 'gradient', true);
    sp2_el = sp_evaluate_element_list (space2, msh_el, 'value', false, 'gradient', true);
    if (nargin == 4)
        for idim = 1:msh.rdim
            x{idim} = reshape (msh_el.geo_map(idim,:,:), msh_el.nqn, msh_el.nel);
        end
        coeffs = coeff (x{:});
    else
        coeffs = ones (msh_el.nqn, msh_el.nel);
    end
    A = A + op_gradu_gradv (sp1_el, sp2_el, msh_el, coeffs);
end
if (nargout == 1)
    varargout{1} = A;
elseif (nargout == 3)
    [rows, cols, vals] = find (A);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
end