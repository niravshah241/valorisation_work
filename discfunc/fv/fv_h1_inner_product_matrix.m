function W = fv_h1_inner_product_matrix(model,model_data)
%function W = fv_h1_inner_product_matrix(model,model_data)
% function computing the h1 inner product matrix for fv-functions on the grid.
%
% By this function, inner products between functions can be
% computed by their DOF vectors: `<U_1, U_2>_{H^1,h} = U_1' W U_2`
%
% return values:
%   W        : inner product matrix

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

grid = model_data.grid;

% currently simple diagonal weighting of element values sufficient
W = sparse(1:grid.nelements,1:grid.nelements,grid.A);
H = grid.nelements;
i=reshape(repmat([1:H],5,1),5*H,1);
j=reshape([1:H;grid.NBI'],5*H,1);
s=100*grid.DC(1,1) * reshape([-ones(1,H);1/4*ones(4,H)],5*H,1);
s(j<1)=0;
j(j<1)=1;
WH = sparse(i,j,s,H,H,5*H);
W = W + WH' * WH;

