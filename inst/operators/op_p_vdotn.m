function varargout = op_p_vdotn(spv, spp, msh)
% ndim = size (gradv, 2);

shpv = reshape (spv.shape_functions, spv.ncomp, msh.nqn, spv.nsh_max, msh.nel);
shpp = reshape (spp.shape_functions, msh.nqn, spp.nsh_max, msh.nel);

ndim = size (shpv, 1);

rows = zeros (msh.nel * spv.nsh_max * spp.nsh_max, 1);
cols = zeros (msh.nel * spv.nsh_max * spp.nsh_max, 1);
values = zeros (msh.nel * spv.nsh_max * spp.nsh_max, 1);

jacdet_weights = msh.jacdet .* msh.quad_weights;

ncounter = 0;
for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))
        p_iel = reshape(shpp(:, :, iel), msh.nqn, spp.nsh_max, 1);
        v_iel = shpv(:, :, :, iel);
        n_iel = reshape (msh.normal(:,:,iel), [ndim, msh.nqn, 1]);
        
        v_cdot_n = reshape(sum(bsxfun(@times, v_iel, n_iel), 1), msh.nqn, spv.nsh_max, 1);
        
        jacdet_iel = reshape (jacdet_weights(:,iel), [msh.nqn,1,1]);

        v_cdot_n_times_j = reshape(bsxfun(@times, v_cdot_n, jacdet_iel), msh.nqn, 1, spv.nsh_max);
        
        tmp = sum(bsxfun(@times, v_cdot_n_times_j, p_iel), 1);
        
        elementary_values = reshape (tmp, spp.nsh_max, spv.nsh_max);
        
        [rows_loc, cols_loc] = ndgrid (spp.connectivity(:,iel), spv.connectivity(:,iel));
        indices = rows_loc & cols_loc;
        rows(ncounter+(1:spp.nsh(iel)*spv.nsh(iel))) = rows_loc(indices);
        cols(ncounter+(1:spp.nsh(iel)*spv.nsh(iel))) = cols_loc(indices);
        values(ncounter+(1:spp.nsh(iel)*spv.nsh(iel))) = elementary_values(indices);
        ncounter = ncounter + spp.nsh(iel)*spv.nsh(iel);
    else
        warning ('geopdes:jacdet_zero_at_quad_node', 'op_p_vdotn: singular map in element number %d', iel)
    end
end

if (nargout == 1 || nargout == 0)
    varargout{1} = sparse (rows, cols, values, spp.ndof, spv.ndof);
elseif (nargout == 3)
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = values;
else
    error ('op_p_vdotn: wrong number of output arguments')
end

end

