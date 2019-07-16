function save_alu3d_hexa(fn,M) 
% function save_alu3d_hexa(fn,M) 
%
% save given ALU3D hexaeder M into alu3D file  
% required fields:
% num_vertices  number of vertices
%     vertices  matrix of vertex coordinates (columnwise)
% num_elements  number of elements
%     elements  matrix of element vertex indices (columnwise)
% num_faces     number of boundary faces
%     faces  matrix of boundary faces (columnwise)

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
  
  fid = fopen(fn,'w');  
  fprintf(fid,'!Hexaeder \n');
  
  % save vertex number
  fprintf(fid,'%d \n',M.num_vertices);
  
  % save vertices
  disp(['writing ',num2str(M.num_vertices),' vertices']);
  fprintf(fid,'%f %f %f \n',M.vertices);
  
  % save element number
  fprintf(fid,'\n %d \n',M.num_elements);
  
  % save elements
  disp(['writing ',num2str(M.num_elements),' elements']);
  fprintf(fid,'%d %d %d %d %d %d %d %d \n',M.elements);
  
  % save boundary segment number
  fprintf(fid,'\n %d \n',M.num_faces);
  
  % save boundary segments
  disp(['writing ',num2str(M.num_faces),' boundary faces']);
  fprintf(fid,'%d %d %d %d %d %d \n',M.faces);
    
  fclose(fid);
  
  
% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
