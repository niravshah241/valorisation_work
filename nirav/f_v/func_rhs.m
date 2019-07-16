function [ res ] = func_rhs( glob,params )
%FUNC_RHS is function for rhs. res column vector.
%glob is global coordinate
res = zeros(params.dimrange,1);

x = glob(1);
y = glob(2);

%% analytical example from paper

% res = [-4*params.kinematic_viscosity(params)*(-1+2*y)*(y^2-6*x*y^2+6*x^2*y^2-y+6*x*y...
%     -6*x^2*y+3*x^2-6*x^3+3*x^4)+1-2*x+4*x^3*y^2*(2*y^2-2*y+1)*(y-1)^2*(-1+2*x)*(x-1)^3 ...
%     4*params.kinematic_viscosity(params)*(-1+2*x)*(x^2-6*x^2*y+6*x^2*y^2-x+6*x*y...
%     -6*x*y^2+3*y^2-6*y^3+3*y^4)+4*x^2*y^3*(-1+2*y)*(y-1)^3*(2*x^2-2*x+1)*(x-1)^2]';

%% Benchmark problem

% res = zeros(params.dimrange,1);

%% Standard

% res = [-4*x^3+2*x+12*params.kinematic_viscosity(params)*y^2-2*params.kinematic_viscosity(params) 0]';
% res = 1e0*[2*params.kinematic_viscosity(params)-1 0]';
res = [0 0]';

end