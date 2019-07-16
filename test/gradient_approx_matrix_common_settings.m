mparams.model_type = 'nonlinear';


params = richards_fv_model(mparams);

params.gridtype = 'rectgrid';
params.tstep = 1;
params.gravity = 0;

params.xrange = [0,1];
params.yrange = [0,1];
params.xnumintervals = 2;
params.ynumintervals = 2;

% set all to dirichlet-boundary by specifying "rectangles", all
% boundary within is set to boundary type by bnd_rect_index
params.bnd_rect_corner1 = [-0.999999, 0; ...
                            0, -0.9999999];
params.bnd_rect_corner2 = [ 2, 0.999999; ...
                            0.999999, 2];

% -1 means dirichlet, -2 means neumann
params.bnd_rect_index =  [ -1, -2];

% name of function in rbmatlab/datafunc/dirichlet_values
params.dirichlet_values_ptr = @dirichlet_values_leftright;

params.c_dir_left    = 1;
params.c_dir_right   = 1;
params.dir_middle    = 0.5;
params.c_dir_correct = 1;
params.decomp_mode   = 0;

params.neumann_values_ptr = @neumann_values_homogeneous;
params.c_neu = 0;

params.rb_init_values                 = @rb_init_values_default;
params.init_values_algorithm          = @fv_init_values;
params.init_values_ptr                = @init_values_homogeneous;
params.implicit_operators_algorithm   = @fv_operators_implicit;
params.operators_diff_implicit        = @fv_operators_diff_implicit_gradient_tensor;
params.operators_neumann_implicit = @fv_operators_zero;
params.fv_impl_conv_weight  = 0;
params.fv_impl_react_weight = 0;
params.fv_impl_diff_weight  = 1.0;
params.get_rb_size          = @get_rb_size;

params.inner_product_matrix_algorithm = @fv_inner_product_matrix;
params.model_type                     = 'nonlinear';
params.diffusivity_ptr        = @diffusivity_homogeneous;
params.diffusivity_tensor_ptr = @diffusivity_tensor_richards;
params.k       = 0.005;

params.geometry_transformation          = 'spline';
params.geometry_spline_type             = 'cubic';
params.hill_height                      = 0.0;
params.geometry_transformation_spline_x = [ 0 0.5 1 ];
params.geometry_transformation_spline_y = [ 0 -0.033 0 ];

params.t     = 0;
params.tstep = 1;
params.dt = 0.1;

params.debug   = true;
params.verbose = 0;

model_data = nonlin_evol_gen_model_data(params);
%| \docupdate 
