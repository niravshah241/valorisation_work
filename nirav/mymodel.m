% grid making and saving
% pdegrid
% save('mygrid','p','e','t')
clc
close all
clear all
%% Parameter-reference generation
disp('Generating reference parameter set')
para1 = 0.7; 
para2 = 0.3; 
parameter_reference_set = [para1 para2];
params.reference_parameter = parameter_reference_set;
disp('Parameter reference set created')


%% Parameter-training generation
disp('Generating training parameter set')
para1 = [0.5 0.8 3]; 
para2 = [0.2 0.4 3]; 
parameter_training_set = gen_parameters( para1,para2);
disp('Parameter training set generation finished')

%% Testing
disp('Generating test parameter set')
para_test_1 = [0.5 0.8 3];
para_test_2 = [0.2 0.4 3];
parameter_test_set = gen_test_parameters(para_test_1,para_test_2);
disp('Test parameter set generated')

%% Online parameter

parameter_online = [0.6 0.4];
params.parameter_online = parameter_online;

%% Actual grid
% ACTUAL GRID DO NOT DELETE

params.mesh_number = 1;
params.gridtype = 'triagrid';
params.grid_initfile = ['mygridnirav', num2str(params.mesh_number), '.mat'];
params.bnd_rect_corner1=[-1,-1;-eps,eps]'; % for analytical
params.bnd_rect_corner2=[eps,1+eps;eps,1-3*10^14*eps]';% for analytical ex.
% params.bnd_rect_corner1=[-1,-1;100,10]'; % for benchmark problem
% params.bnd_rect_corner2=[2,2;100,10-eps]'; % for benchmark problem
% params.bnd_rect_corner1=[-1,-1;1-eps,3*10^14*eps]'; % for standard
% params.bnd_rect_corner2=[eps,1+eps;1+eps,1-eps]'; % for standard
params.bnd_rect_index=[-1,-2];
grid = construct_grid(params);
show_sparsity = false; % Bool variable which plots sparsity pattern of
% assembled matrix is set to true else(i.e. false) the sparsity pattern is not shown
params.show_sparsity = show_sparsity;
paramsP.show_sparsity = show_sparsity;

phase = 'reference';

%ACTUAL GRID OVER

[ grid, params, ref_el_subd_1, ref_el_subd_2, ref_el_subd_3, ...
    ref_el_subd_4, ref_el_subd_5] = grid_dd_func( params, phase );

ref_el_subd{1} = ref_el_subd_1;
ref_el_subd{2} = ref_el_subd_2;
ref_el_subd{3} = ref_el_subd_3;
ref_el_subd{4} = ref_el_subd_4;
ref_el_subd{5} = ref_el_subd_5;

subd_jac{1} = eye(2);
subd_jac{2} = eye(2);
subd_jac{3} = eye(2);
subd_jac{4} = eye(2);
subd_jac{5} = eye(2);

neumann_jac = [1 2 3 4 5];

params.pdeg = 2;
paramsP.pdeg = params.pdeg-1;%taylor hood element
params.dimrange = 2;
paramsP.dimrange = 1;

nrep=[3 6 10 15];

params.ndofs_per_element= nrep(params.pdeg) * params.dimrange;
params.ndofs = params.ndofs_per_element * grid.nelements;
params.dofs = zeros(params.ndofs,1);

paramsP.ndofs_per_element= nrep(paramsP.pdeg) * paramsP.dimrange;
paramsP.ndofs = paramsP.ndofs_per_element * grid.nelements;
paramsP.dofs = zeros(paramsP.ndofs,1);

params.qdeg = params.pdeg;
params.mu = 4;
params.kinematic_viscosity = @(params) params.mu * 1e-6;
mu = params.kinematic_viscosity(params);
c11 = 1e2;% penalty parameter, must be large enough for coercivity

params.rhs_func = @(glob,params,paramsP,grid) [glob(1)  glob(2)]';