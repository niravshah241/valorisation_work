function [ res ] = dirichlet_values(glob,params)
% Function evluating dirichlet boundary velocity value based on pointwise evaluation
% Inputs are global corordinates. Output is column vector
%% For actual case
% if glob(1)<eps
%     res = 100*glob(2)*(1-glob(2));
%     res = [res 0]';
%     
% elseif (glob(2)>(1-eps) || glob(2)<eps)
%     res = [0 0]';
%     
% else
%     res = [0 0]'; %on obstacle
% end

%% for lid driven cavity
% if glob(1) < 0.1 & glob(2) > 1-eps
%     res = [10*glob(1) 0]';
% elseif glob(1) >= 0.1 & glob(1) <= 0.9 & glob(2) > 1-eps
%     res = [1 0]';
% elseif glob(1) > 0.9 & glob(2) > 1-eps
%     res = [10-10*glob(1) 0]';
% else
%     res = [0 0]';
% end

%% for analytical case
if glob(1) < eps
    res = glob(2) * (1 - glob(2));
    res = [res 0]';
else
    res = [0 0]';
end

end