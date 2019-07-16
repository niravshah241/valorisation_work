function [ res ] = pressure_velocity_continuity_integral_21_extension...
    ( params, paramsP, k, grid, qdeg)
%PRESSURE_VELOCITY_CONTINUITY_INTEGRAL_11_EXTENSION Summary of this function goes here
%   Detailed explanation goes here

if nargin == 4
    qdeg = params.qdeg;
end

f = @(lcoord) pressure_velocity_continuity_local_21_extension...
    ( lcoord, params, paramsP, k, grid);

res = -triaquadrature(qdeg,f) * 2 * grid.A(k);

end