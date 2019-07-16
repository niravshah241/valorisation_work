function res=source_integral(params, paramsP, grid, k, qdeg);
%B=[grid.JIT(k,:,1);grid.JIT(k,:,2)];

if nargin == 4
    qdeg = params.qdeg;
end

f=@(lcoord) source(lcoord,params, paramsP, grid, k);
res=triaquadrature(qdeg,f)*2*grid.A(k);%2*grid.A is jacobian determinant
end