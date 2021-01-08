% SP_ERROR_PRESS_DARCY: Evaluate the error in the Nitsche norm scaling as the H^1 error.
%
%   [errl2, varargout] = sp_error_nitsche (space, msh, u, uex, varargin)
%
% INPUT:
%
%   space:       object defining the space of discrete functions (see sp_scalar)
%   msh:         object defining the domain partition and the quadrature rule (see msh_cartesian)
%   u:           vector of dof weights
%   uex:         function handle to evaluate the exact solution
%   varargin:    boundary sides where essetnial b.c. are weakly imposed via Nitsche
%
% OUTPUT:
%
%   errl2:      error in L^2 norm
%   varargout:  error in the Nitsche norm scaling as the L^2 norm
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

function [errl2, varargout] = sp_press_error_darcy (space, msh, u, uex,varargin)

errl2 = 0;

for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (space, msh_col, 'value', true, 'gradient', false);
    
    errl2 = errl2 + (sp_l2_error (sp_col, msh_col, u, uex)).^2;
end

errl2 = sqrt (errl2);

if not(isempty(varargin))
    errl2_bd_square = 0;
    weak_drchlt_sides = varargin{1};
    for iside = weak_drchlt_sides
        msh_side = msh_eval_boundary_side(msh, iside);
        msh_side_from_interior = msh_boundary_side_from_interior (msh, iside);
        sp_bnd = space.constructor (msh_side_from_interior);
        sp_bnd = struct (sp_precompute (sp_bnd, msh_side_from_interior, 'value', true));
        hel = msh_side.charlen(1);
        errl2_bd_square = errl2_bd_square + hel*(sp_l2_error(sp_bnd, msh_side, u, uex)).^2;
    end
end

if nargout > 1
    varargout{1} = sqrt(errl2^2 + errl2_bd_square);
end

end

