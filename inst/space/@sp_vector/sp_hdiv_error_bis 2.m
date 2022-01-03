% SP_HDIV_ERROR: Evaluate the error in H(div) norm.
%
%   [errhdiv, errl2, errhdivs] = sp_hdiv_error (space, msh, u, uex, divuex, varargin);
%
% INPUT:
%
%    space:   object defining the space of discrete functions (see sp_vector)
%    msh:     object defining the domain partition and the quadrature rule (see msh_cartesian)
%    u:       vector of dof weights
%    uex:     function handle to evaluate the exact solution
%    divuex:  function handle to evaluate the div of the exact solution
%    varargin:  boundaries weakly imposed essential b.c.
%
%
% OUTPUT:
%
%     errhdiv:       error in H(div) norm
%     errl2:         error in L^2 norm
%     errhdivs:      error of the div in L^2 norm
%     varargout:     error in Nitsche norm (related to Darcy formulation for weakly imposing essential b.c.)
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

function [errhdiv, errl2, errhdivs, varargout] = sp_hdiv_error_bis(space, msh, u, uex, divuex, varargin)

  errl2 = 0; errhdivs = 0;
  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (space, msh_col, 'value', true, 'divergence', true);
    
    [~, err_l2, errhdivs] = sp_hdiv_error (sp_col, msh_col, u, uex, divuex);
    
    errhdivs = errhdivs + errhdivs.^2;
    errl2 = errl2 + err_l2.^2;
  end
  
  errhdiv = sqrt (errl2 + errhdivs);
  errl2 = sqrt (errl2);
  errhdivs = sqrt (errhdivs);
  
  if not(isempty(varargin))
    weak_drchlt_sides = varargin{1};  
    errl2_dotn_bd_square = 0;
    for iside = weak_drchlt_sides
        msh_side = msh_eval_boundary_side(msh, iside);
        msh_side_from_interior = msh_boundary_side_from_interior(msh, iside);
        sp_bnd = space.constructor (msh_side_from_interior);
        sp_bnd = struct (sp_precompute (sp_bnd, msh_side_from_interior, 'value', true));
        hel = msh_side.charlen(1);
       errl2_dotn_bd_square = errl2_dotn_bd_square + hel*(sp_error_udotn (sp_bnd, msh_side, u, uex)).^2;
    end
  end

  % Nitsche error
if nargout > 3
    varargout{1} = sqrt(errl2^2+errl2_dotn_bd_square);
end

end
