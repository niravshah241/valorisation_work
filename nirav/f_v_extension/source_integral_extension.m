function [ res ] = source_integral_extension( params, paramsP, ...
    grid, k)

f = @(lcoord) source_local_extension(lcoord, params, paramsP, grid, k);
res = triaquadrature(params.qdeg,f) * 2 * grid.A(k);

end