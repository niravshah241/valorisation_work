function Mmerge = merge_alu3d_hexa(M1,M2);
%function Mmerge = merge_alu3d_hexa(M1,M2);
%
% function performing a merging of two alugrid3d hexa meshes. Assuming
% that they are non-overlapping, i.e. no domain-check is performed, but simple
% index-shift

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

  
% Bernard Haasdonk 13.3.2006
  
  clear('Mmerge');
  Mmerge.num_vertices = M1.num_vertices + M2.num_vertices;
  Mmerge.vertices = [ M1.vertices, M2.vertices];
  Mmerge.num_elements = M1.num_elements + M2.num_elements;
  % perform index-shift in element-description for second mesh:
  Mmerge.elements = [ M1.elements, M2.elements + M1.num_vertices];
  Mmerge.num_faces = M1.num_faces + M2.num_faces;
  % perform index-shift in boundary vertex indices for second mesh:
  M2faces_shifted = M2.faces + ...
      repmat([0 0 1 1 1 1]'*M1.num_vertices,1,M2.num_faces); 
  Mmerge.faces = [ M1.faces, M2faces_shifted];
  
  
  
  

% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
