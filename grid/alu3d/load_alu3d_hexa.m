function M = load_alu3d_hexa(fn) 
%function M = load_alu3D_hexa(fn) 
%
% load given ALU3D hexaeder file into mesh-structure M  
% resulting fields:
% num_vertices  number of vertices
%     vertices  matrix of vertex coordinates (columnwise)
% num_elements  number of elements
%     elements  matrix of element vertex indices (columnwise)
% num_faces     number of boundary faces
%     faces     matrix of boundary faces (columnwise) 
%
% note that the numbering of the vertices in elements and faces is 
% C-style, i.e. starting from 0
% for working in matlab, the offset 1 has to be added, e.g.
% z-coordinates of 2nd vertex of element 5 is 
%    vertices(3,M.elements(2,5)+1);
% The boundary faces are coded by 6 integers (for hexamesh)
% \#bndtypeId \#nvertices=4 \#vind1 \#vind2 \#vind3 \#vind4

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
  
  clear('M');
  
  fid = fopen(fn,'r');  
  e = fscanf(fid,'%s',1);
  if ~isequal(e,'!Hexaeder')
    keyboard
    error('wrong file format');
  end;
  
  % load vertex number
  M.num_vertices = fscanf(fid,'%d',1);
  
  % load vertices
  disp(['reading ',num2str(M.num_vertices),' vertices']);
  M.vertices = reshape(fscanf(fid,'%f',M.num_vertices*3),3,M.num_vertices);
  
  % load element number
  M.num_elements = fscanf(fid,'%d',1);
  
  % load elements
  disp(['reading ',num2str(M.num_elements),' elements']);
  M.elements = reshape(fscanf(fid,'%d',...
				M.num_elements * 8),8,M.num_elements);
  
  % load boundary segment number
  M.num_faces = fscanf(fid,'%d',1);
  
  % load boundary segments
  disp(['reading ',num2str(M.num_faces),' boundary faces']);
  M.faces = reshape(fscanf(fid,'%d',M.num_faces*6),6,M.num_faces);
  
%  disp('ignoring occasionally processor-distribution information');
  
  fclose(fid);
  
  
% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
