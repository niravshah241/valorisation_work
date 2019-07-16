function cog = cog_faces_alu3d_hexa(M);
%function cog = cog_faces_alu3d_hexa(M);
%
% function computing the cog of the boundary_faces of a hexaeder alu3d mesh

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
  
  [XX,YY,ZZ] = coord_faces_alu3d_hexa(M);
  X = sum(XX)/4;
  Y = sum(YY)/4;
  Z = sum(ZZ)/4;
  cog = [X;Y;Z];
  
  

% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
