function [ grid, params, el_subd_1, el_subd_2, el_subd_3, ...
    el_subd_4, el_subd_5] = ...
    grid_dd_func( params, phase )

if phase == 'reference'
    mu1 = params.reference_parameter(1);
    mu2 = params.reference_parameter(2);
elseif phase == 'train'
    mu1 = params.parameter_training_set(1);
    mu2 = params.parameter_training_set(2);
elseif phase == 'test'
    mu1 = params.parameter_test_set(1);
    mu2 = params.parameter_test_set(2);
elseif phase == 'online'
    mu1 = params.parameter_online(1);
    mu2 = params.parameter_online(2);
else
    error('This phase is not supported')
end

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
[dl,bt] = decsg(gd,sf,ns);
model = createpde(1);
geometryFromEdges(model,dl);
pdegplot(model,'EdgeLabels','on')
axis tight
generateMesh(model);
pdemesh(model)
axis equal
[p,e,t] = meshToPet(generateMesh(model));
save('mygrid_dd_test','p','e','t');
figure()
el_subd_1 = pdesdt(t,1);
i = el_subd_1;
X = p(1,t(1:3,i));
Y = p(2,t(1:3,i));
plot(X,Y,'*')
title('Subdomain 1')
axis tight
figure()
el_subd_2 = pdesdt(t,2);
i = el_subd_2;
X = p(1,t(1:3,i));
Y = p(2,t(1:3,i));
plot(X,Y,'*')
title('Subdomain 2')
axis tight
figure()
el_subd_3 = pdesdt(t,3);
i = el_subd_3;
X = p(1,t(1:3,i));
Y = p(2,t(1:3,i));
plot(X,Y,'*')
title('Subdomain 3')
axis tight
figure()
el_subd_4 = pdesdt(t,4);
i = el_subd_4;
X = p(1,t(1:3,i));
Y = p(2,t(1:3,i));
plot(X,Y,'*')
title('Subdomain 4')
axis tight
figure()
el_subd_5 = pdesdt(t,5);
i = el_subd_5;
X = p(1,t(1:3,i));
Y = p(2,t(1:3,i));
plot(X,Y,'*')
title('Subdomain 5')
axis tight
figure()
el_subd_6 = pdesdt(t,6);
i = el_subd_6;
X = p(1,t(1:3,i));
Y = p(2,t(1:3,i));
plot(X,Y,'*')
title('Subdomain 6')
axis tight
figure()
el_subd_7 = pdesdt(t,7);
i = el_subd_7;
X = p(1,t(1:3,i));
Y = p(2,t(1:3,i));
plot(X,Y,'*')
title('Subdomain 7')
axis tight
figure()
el_subd_8 = pdesdt(t,8);
i = el_subd_8;
X = p(1,t(1:3,i));
Y = p(2,t(1:3,i));
plot(X,Y,'*')
title('Subdomain 8')
axis tight
figure()
el_subd_9 = pdesdt(t,9);
i = el_subd_9;
X = p(1,t(1:3,i));
Y = p(2,t(1:3,i));
plot(X,Y,'*')
title('Subdomain 9')
axis tight

params.grid_initfile = ['mygridnirav',...
    num2str(params.mesh_number),'_extension.mat'];
params.bnd_rect_corner1=[-1,-1;-eps,eps]';
params.bnd_rect_corner2=[eps,1+eps;eps,1-3*10^14*eps]';
params.bnd_rect_index=[-1,-2];
grid = construct_grid(params);
plot(grid);

end