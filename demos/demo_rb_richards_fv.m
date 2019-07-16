% A small script demonstrating the online phase of the two models from the
% diploma thesis "Reduzierte-Basis Methoden für unges\"{a}ttigte
% Grundwasserstr\"{o}mungen"
%

help demo_rb_richards_fv

disp('==================================================================');
disp(' ');
disp('First example:');
disp('Linear heat equation with parabolic geometry transformation');
disp(' ');
disp('==================================================================');

load richards_fv_detailed_interpol
detailed_data.grid = construct_grid(model);

demo_rb_gui(model,detailed_data,[],plot_params);

disp('Press any key to continue...');
pause

disp('=================================================================');
disp(' ');
disp('Second example:');
disp('Nonlinear Richards equation with affine geometry transformation');
disp(' ');
disp('=================================================================');

load richards_affine_detailed_interpol
detailed_data.grid = construct_grid(model);

demo_rb_gui(model,detailed_data,[],plot_params);

