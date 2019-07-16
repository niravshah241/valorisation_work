function [ res ] = pressure_velocity_continuity_integral( params, paramsP,...
    k, grid, qdeg )
%PRESSURE_VELOCITY_CONTINUITY_LOCAL_INTEGRAL Summary of this function goes
%here: TODO
%   Detailed explanation goes here: TODO

f = @(lcoord) pressure_velocity_continuity_local( lcoord, params, paramsP,...
    k, grid);
res = -triaquadrature(qdeg,f)*2*grid.A(k);

end

