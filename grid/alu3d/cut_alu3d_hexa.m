function Mcut = cut_alu3d_hexa(M,pmin,pmax,bndval)
%function Mcut = cut_alu3d_hexa(M,pmin,pmax,bndval)
%
% function cutting out a "cube" out of given ALU3D grid.
% 1. in a first loop, all elements are identified, the cog of which
%    lies strictly within the specified cube, these are marked for
%    deletion
% 2. the boundary segments of to be deleted elements are marked for deletion
% 3. in a second loop all surviving neighbours of deleted elements
%    a new boundary face is marked for inclusion
% 4. the new element and face list is generated
% 5. all vertices are determined, which are not used
%    in any element or face and marked for deletion
% 6. the new vertex, element and boundary list are condensed

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
  
debug_ = 0;
  
% 1. search all elements with cog inside cube for deletion
%    compute cog of all elements

disp('1. detect elements for deletion');
linvind = M.elements(:)+1;
X = sum(reshape(M.vertices(1,linvind),8,M.num_elements))/8;
Y = sum(reshape(M.vertices(2,linvind),8,M.num_elements))/8;
Z = sum(reshape(M.vertices(3,linvind),8,M.num_elements))/8;
cog = [X;Y;Z];

%    must be positive for points within the cube
pmindiff = min(cog-repmat(pmin(:),1,M.num_elements));
pmaxdiff = min(repmat(pmax(:),1,M.num_elements)-cog);
mindiff = min([pmindiff; pmaxdiff]);
el_to_be_del_flag = (mindiff>0);
el_to_be_del_ind = find(mindiff>0);

if debug_
  disp('halt after el_to_be_del detection')
  keyboard;
end;

% 2. mark boundary segments of to be deleted elements for deletion/
%    generate new boundary

% generate list of all faces of to be deleted elements
% in ALUgrid notation 6 faces of 8 vertices:
% modified indices  
%1 2 3 4 5 6 7 8 
%  =>
% 1 4 8 5 
% 2 6 7 3 
% 1 5 6 2 
% 4 3 7 8 
% 1 2 3 4 
% 6 5 8 7 

disp('2. make boundary list of to be deleted elements');
fv_el_ind = [ 1 4 8 5 ...
	  2 6 7 3 ...
	  1 5 6 2 ...
	  4 3 7 8 ...
	  1 2 3 4 ...
	  6 5 8 7];

fv_ind = M.elements(fv_el_ind,el_to_be_del_ind)+1;
fv_ind = reshape(fv_ind,4,6*length(el_to_be_del_ind));
% fvind now is the list of columnwise indices of candidate boundary faces

if debug_
  disp('halt after setup of list of candidate boundary faces')
  keyboard;
end;

% generate description of all faces and process 
%   1-neigh to be del            => do nothing
%   2-neigh not to be del        => insert new boundary
%   3-neigh is existing boundary => delete boundary from existing list
%   0-not qualified neighbour: should not happen

% initially nothing is known about faces
fv_ind_descr = zeros(1,size(fv_ind,2));
% detection of case 1: double appearance of the face in the list:
%keyboard;
disp('2. a) search boundaries to elements which are also to be deleted');
for facenr = 1:length(fv_ind_descr)
  if mod(facenr,100)==0
    disp([num2str(facenr),'/',num2str(length(fv_ind_descr))]);
  end;
  % search for faces with at least one identical vertex nr
  cand_flag1 = max(fv_ind==fv_ind(1,facenr));
  cand_flag2 = max(fv_ind==fv_ind(2,facenr));
  cand_flag3 = max(fv_ind==fv_ind(3,facenr));
  cand_flag4 = max(fv_ind==fv_ind(4,facenr));
  cand_flag = cand_flag1 | cand_flag2 | cand_flag3 | cand_flag4;
  ind_cand = find(cand_flag==1);
  matrix = fv_ind(:,ind_cand);
  ind = face_in_matrix(fv_ind(:,facenr),matrix);
  if (length(ind)>1)  % at least once as identical
    fv_ind_descr(facenr) = 1; 
  end;
end;

if debug_
  disp('halt after detection of case 1 of boundary faces')
  keyboard;
end;

% detection of case 2: 
%   generate element list of elements with any vertex appearing in
%   tobedeleted elements.
disp('2. b) search boundaries to elements which are not be deleted');

vind = zeros(size(M.elements));
el_to_be_del_vind = M.elements(:,el_to_be_del_ind)+1;
j=0;
for i = el_to_be_del_vind(:)'
  j=j+1;
  if mod(j,100)==0
    disp([num2str(j),'/',num2str(length(el_to_be_del_vind))]);
  end;
  fi = find(M.elements == (i-1));
  vind(fi) = 1;
end;
cand_el_flag = max(vind); 

%   make sure, that all listed elements are surviving
cand_el_flag(el_to_be_del_ind) = 0;
cand_el_ind = find(cand_el_flag==1);

if debug_
  disp('halt after detection of surviving neighbour elements')
  keyboard;
end;

%   generate all boundary faces of remaining elements
fv_boundary_ind = M.elements(fv_el_ind,cand_el_ind)+1;
fv_boundary_ind = reshape(fv_boundary_ind,4,6*length(cand_el_ind));

%   search all questionable faces in the generated list
face_to_be_insert_flag = zeros(1,size(fv_boundary_ind,2));
for facenr = 1:length(fv_ind_descr)
  if mod(facenr,100)==0
    disp([num2str(facenr),'/',num2str(length(fv_ind_descr))]);
  end;
  % search for faces with at least one identical vertex nr
  cand_flag1 = max(fv_boundary_ind==fv_ind(1,facenr));
  cand_flag2 = max(fv_boundary_ind==fv_ind(2,facenr));
  cand_flag3 = max(fv_boundary_ind==fv_ind(3,facenr));
  cand_flag4 = max(fv_boundary_ind==fv_ind(4,facenr));
  cand_flag = cand_flag1 | cand_flag2 | cand_flag3 | cand_flag4;
  ind_cand = find(cand_flag==1);
  matrix = fv_boundary_ind(:,ind_cand);
  ind = face_in_matrix(fv_ind(:,facenr),matrix);
  if (length(ind)>0)  %
    fv_ind_descr(facenr) = 2; 
  end;
  face_to_be_insert_flag(ind_cand(ind)) = 1;
end;
face_to_be_insert_ind = find(face_to_be_insert_flag==1);

if debug_
  disp('halt after detection of case 2 of boundary faces')
  keyboard;
end;

% detection of case 3: 
disp('2. c) search boundaries to former domain boundary');
face_to_be_del_flag = zeros(1,size(M.faces,2));
for facenr = 1:length(fv_ind_descr)
  if mod(facenr,100)==0
    disp([num2str(facenr),'/',num2str(length(fv_ind_descr))]);
  end;
  % search for faces with at least one identical vertex nr
  cand_flag1 = max(M.faces(3:6,:)== (fv_ind(1,facenr)-1));
  cand_flag2 = max(M.faces(3:6,:)== (fv_ind(2,facenr)-1));
  cand_flag3 = max(M.faces(3:6,:)== (fv_ind(3,facenr)-1));
  cand_flag4 = max(M.faces(3:6,:)== (fv_ind(4,facenr)-1));
  cand_flag = cand_flag1 | cand_flag2 | cand_flag3 | cand_flag4;
  ind_cand = find(cand_flag==1);
  matrix = M.faces(3:6,ind_cand)+1;
  ind = face_in_matrix(fv_ind(:,facenr),matrix);
  if (length(ind)>0)  %
    fv_ind_descr(facenr) = 3; 
  end;
  face_to_be_del_flag(ind_cand(ind)) = 1;
end;
face_to_be_del_ind = find(face_to_be_del_flag==1);

if debug_
  disp('halt after detection of case 3 of boundary faces')
  keyboard;
end;

% make sure that no face is left out of assignment 1-3
if ~isempty(find(fv_ind_descr==0))
  disp('not all questionable boundaries are qualified!');
  keyboard;
end;

disp('3. c) generate new element, faces and vertices');
% 4. generate new element and face list
%    new elements:
Mcut.elements = M.elements(:,find(el_to_be_del_flag==0));
Mcut.num_elements = size(Mcut.elements,2);
%    new faces:
facesnew = [bndval * ones(1,length(face_to_be_insert_ind)); ...
	    4      * ones(1,length(face_to_be_insert_ind)); ...
	    fv_boundary_ind(:,face_to_be_insert_ind)-1 ];
Mcut.faces = [M.faces(:,find(face_to_be_del_flag==0)), ...
	      facesnew];
Mcut.num_faces = size(Mcut.faces,2);

if debug_
  disp('halt after setup of non-condensed result')
  keyboard;
end;

% 5. all vertices are determined, which are to be reused
new2old_vflag = zeros(1,M.num_vertices);
new2old_vflag(Mcut.elements+1) = 1;
new2old_vflag(Mcut.faces(3:6,:)+1) = 1;
new2old_vind = find(new2old_vflag==1);
Mcut.num_vertices = length(new2old_vind);
old2new_vind = zeros(1,M.num_vertices);
for i=1:length(new2old_vind);
  old2new_vind(new2old_vind(i)) = i;
end;

% 6. the new vertex, element and boundary list are condensed
Mcut.vertices = M.vertices(:,new2old_vind);
Mcut.elements = reshape(old2new_vind(Mcut.elements(:)+1),...
			8,Mcut.num_elements)-1;
Mcf = Mcut.faces(3:6,:);
Mcut.faces(3:6,:) = reshape(old2new_vind(Mcf(:)+1),4,Mcut.num_faces)-1;

if debug_
  disp('halt before finishing routine');
  keyboard;
end;


% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
