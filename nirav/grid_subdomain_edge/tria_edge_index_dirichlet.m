function [ tria_index,local_edge_index,a_index, b_index, local_vertex_index] = tria_edge_index_dirichlet( grid )
%TRIA_EDGE_INDEX_DIRICHLET Summary of this function goes here:TODO
%   Detailed explanation goes here:TODO

[tria_index,local_vertex_index]=find(grid.NBI==-1);

tria_index=tria_index'; %transposed only to print better on screen
local_vertex_index=local_vertex_index';

for i=1:1:length(tria_index)
    local_edge_index(i)=grid.INB(tria_index(i),local_vertex_index(i));
    a_index(i)=grid.VI(tria_index(i),local_vertex_index(i));
end

b_index=a_index+1;

figure()
plot(grid.X(a_index),grid.Y(a_index),'.');
axis([min(grid.X)-0.1 max(grid.X)+0.1 min(grid.Y)-0.1 max(grid.Y)+0.1]);
title('Dirichlet boundary');

end