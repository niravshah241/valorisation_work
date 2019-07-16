function [C, G] = aff_trafo_loc2glob(x0, y0)
%function [C, G] = aff_trafo_loc2glob(x0, y0)
% function giving the coefficients for the affine transformation
% from reference/local triangle to the original/global one.
%
% In detail, this implements a transformation of type
% ``T^{-1}_{i,aff}(x;\mu) = C^{-1}_{i,aff}(\mu) + \sum_{j=1,2} G^{k,-1}_{ij}(\mu) x_j \qquad      i=1,2``
%
% The transformation is calculated for the standard triangle:
%
% @verbatim
% triangle:   (0,1)  |\                                   /| (x0(3),y0(3))
%                    | \        ---T--->                 / |
%             (0,0)  |__\ (1,0)           (x0(1),y0(1)) /__| (x0(2),y0(2))
% @endverbatim
%
% Parameters:
%  x0: vector of size '3 x 1' holding x values of the original/global triangle
%  y0: vector of size '3 x 1' holding y values of the original/global triangle
%
% Return values:
%  C: matrix of size '2 x 1' with entries 'C=[c1; c2]'
%  G: matrix of size '2 x 2' with entries 'G=[g11, g12; g21, g22]'
%
% See also aff_trafo_glob2loc() which gives the transformation in the other
% direction (global to local)

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


% Oliver Zeeb, 01.02.11


%standard triangle
x1 = 0;
y1 = 0;
x2 = 1;
y2 = 0;
x3 = 0;
y3 = 1;

    

%matrices for transformation:
B_aff =  [1, 0, x1, y1, 0, 0; ...
          0, 1, 0, 0, x1, y1; ...
          1, 0, x2, y2, 0, 0; ...
          0, 1, 0, 0, x2, y2; ...
          1, 0, x3, y3, 0, 0; ...
          0, 1, 0, 0, x3, y3];
     
V_aff = [x0(1); y0(1); x0(2); y0(2); x0(3); y0(3)];

coef_vec = B_aff \ V_aff;

C = [coef_vec(1); coef_vec(2)];
G = [coef_vec(3), coef_vec(4); coef_vec(5), coef_vec(6)];
% c1  = coef_vec(1);
% c2  = coef_vec(2);
% g11 = coef_vec(3);
% g12 = coef_vec(4);
% g21 = coef_vec(5);
% g22 = coef_vec(6);
