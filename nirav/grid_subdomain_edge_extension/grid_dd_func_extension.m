function [ grid, params, el_subd] = ...
    grid_dd_func_extension(params)

mu1 = params.reference_parameter(1);
mu2 = params.reference_parameter(2);

triangle1x = [0 0 0.1];
triangle1y = [0 1 1];

triangle2x = [0 0.1 0.1];
triangle2y = [0 0 1];

triangle3x = [0.1 0.3 mu1];
triangle3y = [0 0 mu2];

triangle4x = [0.1 0.1 mu1];
triangle4y = [0 1 mu2];

triangle5x = [0.1 0.9 mu1];
triangle5y = [1 1 mu2];

triangle6x = [0.9 0.9 mu1];
triangle6y = [0 1 mu2];

triangle7x = [0.9 0.7 mu1];
triangle7y = [0 0 mu2];

triangle8x = [0.9 0.9 1];
triangle8y = [0 1 0];

triangle9x = [0.9 1 1];
triangle9y = [1 1 0];

T1 = [2 3 triangle1x triangle1y]';
T2 = [2 3 triangle2x triangle2y]';
T3 = [2 3 triangle3x triangle3y]';
T4 = [2 3 triangle4x triangle4y]';
T5 = [2 3 triangle5x triangle5y]';
T6 = [2 3 triangle6x triangle6y]';
T7 = [2 3 triangle7x triangle7y]';
T8 = [2 3 triangle8x triangle8y]';
T9 = [2 3 triangle9x triangle9y]';

gd = [T1,T2,T3,T4,T5,T6,T7,T8,T9];
ns = char('T1','T2','T3','T4','T5','T6','T7','T8','T9');
ns = ns';
sf = '(T1+T2+T3+T4+T5+T6+T7+T8+T9)';
[dl, bt] = decsg(gd,sf,ns);
model = createpde(1);
geometryFromEdges(model,dl);
%figure()
%pdegplot(model,'EdgeLabels','on')
%axis tight
%figure()
%generateMesh(model);
%axis tight
%figure()
%pdemesh(model)
%axis tight
%axis equal
[p,e,t] = meshToPet(generateMesh(model));
save(strcat('mygridnirav',...
    num2str(params.mesh_number),'_extension'),'p','e','t');
el_subd_1 = pdesdt(t,1);
el_subd_2 = pdesdt(t,2);
el_subd_3 = pdesdt(t,3);
el_subd_4 = pdesdt(t,4);
el_subd_5 = pdesdt(t,5);
el_subd_6 = pdesdt(t,6);
el_subd_7 = pdesdt(t,7);
el_subd_8 = pdesdt(t,8);
el_subd_9 = pdesdt(t,9);

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
    figure()
    i = el_subd_3;
    X = p(1,t(1:3,i));
    Y = p(2,t(1:3,i));
    plot(X,Y,'*')
    title('Subdomain 3')
    axis tight
    axis equal
    figure()
    i = el_subd_4;
    X = p(1,t(1:3,i));
    Y = p(2,t(1:3,i));
    plot(X,Y,'*')
    title('Subdomain 4')
    axis tight
    axis equal
    figure()
    i = el_subd_5;
    X = p(1,t(1:3,i));
    Y = p(2,t(1:3,i));
    plot(X,Y,'*')
    title('Subdomain 5')
    axis tight
    axis equal
    figure()
    i = el_subd_6;
    X = p(1,t(1:3,i));
    Y = p(2,t(1:3,i));
    plot(X,Y,'*')
    title('Subdomain 6')
    axis tight
    axis equal
    figure()
    i = el_subd_7;
    X = p(1,t(1:3,i));
    Y = p(2,t(1:3,i));
    plot(X,Y,'*')
    title('Subdomain 7')
    axis tight
    axis equal
    figure()
    i = el_subd_8;
    X = p(1,t(1:3,i));
    Y = p(2,t(1:3,i));
    plot(X,Y,'*')
    title('Subdomain 8')
    axis tight
    axis equal
    figure()
    i = el_subd_9;
    X = p(1,t(1:3,i));
    Y = p(2,t(1:3,i));
    plot(X,Y,'*')
    title('Subdomain 9')
    axis tight
    axis equal
end

params.grid_initfile = ['mygridnirav',...
    num2str(params.mesh_number),'_extension.mat'];
params.bnd_rect_corner1=[-1,-1;1-eps,3*10^14*eps]'; % for standard
params.bnd_rect_corner2=[eps,1+eps;1+eps,1-eps]'; % for standard
params.bnd_rect_index=[-1,-2];
grid = construct_grid(params);
if params.grid_show == true
    figure()
    plot(grid);
end

el_subd = { el_subd_1, el_subd_2, el_subd_3, el_subd_4, ...
    el_subd_5, el_subd_6, el_subd_7, el_subd_8, el_subd_9};

a1 = size(el_subd);

num_elements = 0;
for i = 1:1:numel(el_subd)
    num_elements = num_elements + numel(el_subd{i});
end

assert(num_elements == grid.nelements,...
    'Number of elements in grid not same as total elements in subdomains.');

end