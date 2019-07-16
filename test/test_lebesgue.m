function OK=test_lebesgue
% this is a script showing the ei_detailed construction for a function which
% empirically interpolated turns out to have the worst possible Lebesgue
% constant `\Lambda = \max_{x} \sum_{m=1}^M |\xi_m(x)| = 2^M - 1`
%

% This program is open source.  For license terms, see the COPYING file.
%
% --------------------------------------------------------------------
% ATTRIBUTION NOTICE:
% This product includes software developed for the RBmatlab project at
% (C) Universities of Stuttgart and Münster, Germany.
%
% RBmatlab is a MATLAB software package for model reduction with an
% emphasis on Reduced Basis Methods. The project is maintained by
% M. Dihlmann, M. Drohmann, B. Haasdonk, M. Ohlberger and M. Schaefer.
% For Online Documentation and Download we refer to www.morepas.org.
% --------------------------------------------------------------------


maxfuncs = 10;

model.gridtype = 'rectgrid';
model.xrange = [ 0, 1];
model.yrange = [ 0, 1];
model.xnumintervals = 100;
model.ynumintervals = 1;
model.verbose       = 1;

model_data.grid = construct_grid(model);
model_data.W    = fv_inner_product_matrix(model, model_data);

model.linfty_error_sequence_algorithm = @fv_linfty_error;

U = zeros(model_data.grid.nelements, maxfuncs);

U(1,:) = 1-100000*eps;

for i=1:maxfuncs
  U(1+2*i,i) = 1;
  U(3+2*i:2:100,i) = -1;
end

params.ei_Mmax = 3;
params.ei_stop_on_Mmax = 3;
params.ei_target_error = 'linfty-interpol';
params.compute_lebesgue = true;

LU = U;
save(fullfile(rbmatlabtemp,'LUtmp'), 'LU');

LU_fnames = {'LUtmp'};

model.get_inner_product_matrix = @(model_data) model_data.W;
model.nt = maxfuncs-1;

detailed_data = ei_detailed(model, model_data, LU_fnames, params);

OK = (detailed_data.ei_info{1}.lebesgue == detailed_data.ei_info{1}.max_lebesgue);

end
