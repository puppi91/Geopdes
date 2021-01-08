% SP_ERROR_PRESS_DARCY: Evaluate the error for the pressure for Darcy
%
%   [err, errh1s, err_jump, err_bd] = sp_error_press_darcy (space, msh, u, uex, graduex, varargin)
%
% INPUT:
%
%   space:       object defining the space of discrete functions (see sp_scalar)
%   msh:         object defining the domain partition and the quadrature rule (see msh_cartesian)
%   u:           vector of dof weights
%   uex:         function handle to evaluate the exact solution
%   graduex:     function handle to evaluate the gradient of the exact solution
%   solver:      either 1 or 2
%   bnd_sides:   boundary sides where essetnial b.c. are weakly imposed via Nitsche

%
% OUTPUT:
%
%   err:      error in suitable norm
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

function [err, err_h1s, err_l2_jump, err_l2_bd, varargout] = sp_error_press_darcy (space, msh, u, uex, graduex, solver, bnd_sides)

p = space.degree(1);
% element size h
msh_elem = msh_evaluate_element_list(msh, 1);
h_el = msh_elem.element_size(1);

% H^1-seminorm contribution
msh_el = msh_precompute (msh);
sp_el  = sp_precompute (space, msh_el, 'gradient', true);
[~,~,errh1s_temp] = sp_h1_error (sp_el, msh_el, u, uex, graduex);
err_h1s = (errh1s_temp).^2;


% L^2-norm of the jumps through the internal faces contribution
err_l2_jump = 0;
for ii = 1:msh.ndim
    for jj = 1: numel(msh.breaks{ii})-2
        [~, msh_left, msh_right] =  msh_on_internal_face (msh, ii, jj+1);
        sp_left = space.constructor (msh_left);
        sp_left = struct (sp_precompute (sp_left, msh_left));
        
        sp_right = space.constructor (msh_right);
        sp_right = struct (sp_precompute (sp_right, msh_right));
        
        err_l2_jump = err_l2_jump + ((h_el).^(-1))*(sp_l2_jump_error(sp_left,...
            msh_left, sp_right, msh_right, u, uex)).^2;        
    end
end

err_l2_bd = 0;
if solver==2
    % L^2-norm on the weak Neumann boundary
    for iside = bnd_sides
        msh_side = msh_eval_boundary_side(msh, iside);
        msh_side_from_interior = msh_boundary_side_from_interior (msh, iside);
        sp_bnd = space.constructor (msh_side_from_interior);
        sp_bnd = struct (sp_precompute (sp_bnd, msh_side_from_interior, 'value', true));
        err_l2_bd = err_l2_bd + ((h_el).^(-1))*(sp_l2_error(sp_bnd, msh_side, u, uex)).^2;
    end
end

err = sqrt(err_h1s+err_l2_jump+err_l2_bd);
err_h1s = sqrt(err_h1s);
err_l2_jump = sqrt(err_l2_jump);
err_l2_bd = sqrt(err_l2_bd);

end

