function cog = cog_elements_alu3d_hexa(M);
%function cog = cog_elements_alu3d_hexa(M);
%
% function computing the cog of the elements of a hexaeder alu3d mesh

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

  
% Bernard Haasdonk 16.3.2006

linvind = M.elements(:)+1; % matlab/C-offset
X = sum(reshape(M.vertices(1,linvind),8,M.num_elements))/8;
Y = sum(reshape(M.vertices(2,linvind),8,M.num_elements))/8;
Z = sum(reshape(M.vertices(3,linvind),8,M.num_elements))/8;
cog = [X;Y;Z];


% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
