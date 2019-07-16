function [XX,YY,ZZ] = coord_faces_alu3d_hexa(M)
%function [XX,YY,ZZ] = coord_faces_alu3d_hexa(M)
%
% function collecting the coordinates of the face patches
% the columns of XX,YY,ZZ represent the x,y,z coordinates of all 4 face points
% result is therefore 3 matrices of size 4 x M.num_faces
% the resulting matrices can directly be used in 'patch' for example.

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
  
  faces_el_ind = M.faces(3:6,:);
  linvind = faces_el_ind(:)+1; % matlab/C-offset
  XX = reshape(M.vertices(1,linvind),4,M.num_faces);
  YY = reshape(M.vertices(2,linvind),4,M.num_faces);
  ZZ = reshape(M.vertices(3,linvind),4,M.num_faces);
  
  
  
% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
