function K = fv_inner_product(model,model_data,U1,U2)
%function K = fv_inner_product(model,[model_data],U1,U2)
%
% function computing the l2 inner product between all pairs of 
% fv-functions in U1,U2. Result is a matrix K of size size(U1,2) x
% size(U2,2) with the inner products.
% Correct Omega-integrals are computed by respecting the cell-areas
% defined in grid,params. Actually, params is currently
% superfluous, but kept for consistency of command line arguments.
% perhaps later params may indicate polynomial degree, etc.

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


% Bernard Haasdonk 20.7.2006

% currently simple diagonal weighting of element values sufficient
if(isempty(model_data))
  model_data = model.gen_model_data(model);
end
W = model.get_inner_product_matrix(model_data);
K = U1' * W * U2;

%| \docupdate 
