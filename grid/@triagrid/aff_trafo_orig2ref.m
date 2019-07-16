function [C, G] = aff_trafo_orig2ref(x0, y0, varargin)
%function [C, G] = aff_trafo_orig2ref(x0, y0, varargin)
% function giving the coefficients for the affine transformation
% from original triangle to the reference one,
%
% ``T^{-1}_{i,aff}(x;\mu) = C^{-1}_{i,aff}(\mu) + \sum_{j=1,2} G^{k,-1}_{ij}(\mu) x_j \qquad      i=1,2``
%
% @verbatim
% triangle:                     /| (x0(3),y0(3))              (0,1)  |\
%                              / |                ---T--->           | \
%               (x0(1),y0(1)) /__| (x0(2),y0(2))              (0,0)  |__\ (1,0)
% @endverbatim
%
% @deprecated I guess, that this function is deprecated and
% aff_trafo_glob2loc() should be used instead...
%
% function giving c1, c2, g11, g12, g21, g22
% so: C=[c1; c2] and G=[g11, g12; g21, g22]

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

x1 = 0;
y1 = 0;
x2 = 1;
y2 = 0;
x3 = 0;
y3 = 1;

ref_coord = [x1; y1; x2; y2; x3; y3];

B_aff =  [1, 0, x0(1), y0(1), 0, 0; ...
          0, 1, 0, 0, x0(1), y0(1); ...
          1, 0, x0(2), y0(2), 0, 0; ...
          0, 1, 0, 0, x0(2), y0(2); ...
          1, 0, x0(3), y0(3), 0, 0; ...
          0, 1, 0, 0, x0(3), y0(3)];
      
coef_vec = B_aff \ ref_coord;

C = [coef_vec(1); coef_vec(2)];
G = [coef_vec(3), coef_vec(4); coef_vec(5), coef_vec(6)];
% c1  = coef_vec(1);
% c2  = coef_vec(2);
% g11 = coef_vec(3);
% g12 = coef_vec(4);
% g21 = coef_vec(5);
% g22 = coef_vec(6);
