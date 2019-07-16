function [S1_without_S2, i1] = vectorset_difference(S1,S2)
%function [S1_without_S2, i1] = vectorset_difference(S1,S2)
%
% function determining the set-theoretic difference 'S1\\S2' of the
% columns in S1 and S2. i1 is a set of indices, which represent the
% columns of S1, i.e. 'S1(:,i1) = S1_without_S2'

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


% Bernard Haasdonk 5.6.2007

% search list of correspondences

[i,j] = find_corresp(S1,S2);

i_mask = zeros(1,size(S1,2));
i_mask(i) = 1;

i1 = find(i_mask==0); 
S1_without_S2 = S1(:,i1);






% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
