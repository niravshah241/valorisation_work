function [ res ] = neumann_values(glob,params,paramsP)
%NEUMANN_VALUES Summary of this function goes here
%   Detailed explanation goes here

%% analytical example
if glob(1) > 1 - eps
    res = [0 0]';
else
    disp('Not on Neumann boundary')
end

%% benchmark problem lid driven cavity
% error('No Neumann boundary should exist here');

%% For actual case
% res = zeros(2,1);
% if glob(1)>(1-eps)
%     res = [0 0]';
% else
%     disp('Not on neumann boundary');
%     %error('Not on neumann boundary');
% end
end