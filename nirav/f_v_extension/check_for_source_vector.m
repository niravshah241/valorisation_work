%TODO : affine transformation of source vector works only for constant
%rhs_func and not when rhs_func is dependent on x,y.

source_vector = source_assembly_extension( params, paramsP, ...
    transformed_grid);
source_vector_real = source_assembly( params, paramsP, transformed_grid);
error_assembly_source_vector = ...
    max(abs(source_vector_real - source_vector));

tic();
source_vector = source_assembly_extension( params, paramsP, grid);
t_end = toc();
disp(['Time for source_assembly_extension : ' num2str(t_end)]);
tic();
source_vector_online = source_assembly_extension_online...
    ( source_vector, params, grid, para_mapping, el_subd);
t_end = toc();
disp(['Time for source_assembly_extension_online : ' num2str(t_end)]);
source_vector_real = source_assembly_extension...
    ( params, paramsP, transformed_grid);
error_affine_source_vector = ...
    max(abs(source_vector_real - source_vector_online));