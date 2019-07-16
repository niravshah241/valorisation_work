function [res] = gen_parameters( para1,para2)
%GEN_TRAIN_PARAMETERS Summary of this function goes here
% %   Detailed explanation goes here

para1 = linspace(para1(1),para1(2),para1(3));
para2 = linspace(para2(1),para2(2),para2(3));
[para1,para2] = meshgrid(para1,para2);
para1 = reshape(para1,[],1);
para2 = reshape(para2,[],1);

res = [para1,para2];

end