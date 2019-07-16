function micro2macro = micro2macro_map(microgrid, macrogrid)
%function micro2macro = micro2macro_map(model, grid)
%
% function defining a vector micro2macro containing the information
% which triangle of the microgrid lies in which triangle of the
% macrogrid, defined in the model
% micro2macro(5) = 7 means that micro-triangle nr 5 
% lies in macro-triangle nr 7
%
% microgrid and macrogrid must be triagrid
%
% Oliver Zeeb, 01.02.11

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


micro2macro = zeros(microgrid.nelements,1);

pmacro_x = macrogrid.X;
pmacro_y = macrogrid.Y;
tmacro = macrogrid.VI;
nr_macro_tri =macrogrid.nelements; %nr of macro-triangles

pmicro_x = microgrid.X;
pmicro_y = microgrid.Y;

%dummy matrices for the affine transfomation from original macro triangle
%to reference triangle
C=zeros(2,nr_macro_tri);
G=zeros(2,2,nr_macro_tri);

% get all the transformations of the macrotriangles to the
% reference-triangle:
for k=1:nr_macro_tri
    tria_pts_x = pmacro_x(tmacro(k,:),:);
    tria_pts_y = pmacro_y(tmacro(k,:),:);
    [C(:,k), G(:,:,k)] = aff_trafo_orig2ref(tria_pts_x, tria_pts_y);
end

%check which point is in which macrotriangle
for macro_element_id = 1:nr_macro_tri
    pts_in_mac_tri = zeros(1,microgrid.nvertices);
    C_big = repmat(C(:,macro_element_id),1,microgrid.nvertices);
    pts_ref_dom = C_big + G(:,:,macro_element_id)*[pmicro_x'; pmicro_y']; %transform all points
    %check, which of the transformed points are in the reference triangle:
    % CAREFUL!!! See the "eps"! The comparison, wheter a value ist bigger 0
    % or smaller 1 is a bit sloppy...
    i=(pts_ref_dom(1,:)>=0-10*eps & pts_ref_dom(2,:)>=0-10*eps & pts_ref_dom(1,:)+pts_ref_dom(2,:)<=1+10*eps);
    %original (mathematically correct, but unfortunatelly not working
    %correctly...):
    %i=(pts_ref_dom(1,:)>=0 & pts_ref_dom(2,:)>=0 & pts_ref_dom(1,:)+pts_ref_dom(2,:)<=1);
    pts_in_mac_tri(i)=1; %index of the points in macro-triangle k
    
    %check, for which microtriangle all the vertices are in the macro-tringle:
    % if all 3 vertices of the transformed mirotriangle are in the
    % standard-triangle, then this microtriangle lies in the macrotriangle
    bool_mat = zeros(size(microgrid.VI));
    bool_mat(:,:) = pts_in_mac_tri(microgrid.VI(:,:));
    %bool_mat is a matrix of triangles with entry = 1 if the corresponding point
    %is in the reference triangle, so the original point is in the
    %macro-triangle and entry = 0 if the correspondong point is not in
    %the macro tiangle
    mic_in_mac_triangle=(bool_mat(:,1)==1 & bool_mat(:,2)==1 & bool_mat(:,3)==1);
    % entry of mic_in_mac_triangle is 1, if all three vertices are in the standard triangle
    % else it is 0 
    micro2macro = micro2macro + macro_element_id.* mic_in_mac_triangle; 
end

