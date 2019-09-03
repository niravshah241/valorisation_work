function [ grid, params, el_subd ] = ...
    grid_dd_func_extension_lid_driven_cavity( params)
%GRID_DD_FUNC_EXTEN Summary of this function goes here
%   Detailed explanation goes here

triangle1x = [0 1 0];
triangle1y = [0 0 1];

triangle2x = [1 0 1];
triangle2y = [0 1 1];

T1 = [2 3 triangle1x triangle1y]';
T2 = [2 3 triangle2x triangle2y]';

gd = [T1,T2];
ns = char('T1','T2');
ns = ns';
sf = '(T1+T2)';
[dl, bt] = decsg(gd,sf,ns);
model = createpde(1);
geometryFromEdges(model,dl);
[p,e,t] = meshToPet(generateMesh(model));
save(strcat('mygridnirav',...
    num2str(params.mesh_number),'_extension_lid_driven_cavity'),'p','e','t');
el_subd_1 = pdesdt(t,1);
el_subd_2 = pdesdt(t,2);
if params.grid_show == true
    figure()
    i = el_subd_1;
    X = p(1,t(1:3,i));
    Y = p(2,t(1:3,i));
    plot(X,Y,'*')
    title('Subdomain 1')
    axis tight
    axis equal
    figure()
    i = el_subd_2;
    X = p(1,t(1:3,i));
    Y = p(2,t(1:3,i));
    plot(X,Y,'*')
    title('Subdomain 2')
    axis tight
    axis equal
end

params.grid_initfile = ['mygridnirav',...
    num2str(params.mesh_number),'_extension_lid_driven_cavity.mat'];
params.bnd_rect_corner1=[-1,-1]'; % for lid driven cavty
params.bnd_rect_corner2=[1+eps,1+eps]'; % for lid driven cavty
params.bnd_rect_index=[-1];
grid = construct_grid(params);
if params.grid_show == true
    figure()
    plot(grid);
end

el_subd = { el_subd_1, el_subd_2};

a1 = size(el_subd);

num_elements = 0;
for i = 1:1:numel(el_subd)
    num_elements = num_elements + numel(el_subd{i});
end

assert(num_elements == grid.nelements,...
    'Number of elements in grid not same as total elements in subdomains.');

end