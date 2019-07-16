function p = plot_bnd_alu3d_hexa(M, params)
%function p = plot_bnd_alu3d_hexa(M)
%
% plot boundary of alu3d-mesh 
% params can contain fields with additional parameters
%      shrink_factor : boundary patches are shrinked with this factor
%                     (default is 1 == no shrinking, <1 is shrinking)
%      face_value    : faces are plotted with a color according to this
%                      vector. If not given, the
%                      boundary-type-assignments are used.

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

    
% Bernard Haasdonk 15.3.2006
  
  if ~isfield(params,'shrink_factor')
    params.shrink_factor = 1;
  end;

  if ~isfield(params,'face_value')
    params.face_value = M.faces(1,:);
  end;
  
  % compute cog of faces
  cog = cog_faces_alu3d_hexa(M);
  
  % compute vertices of faces
  [XX,YY,ZZ] = coord_faces_alu3d_hexa(M);
  
  % perform shrinking 
  cogX = repmat(cog(1,:),4,1);
  XX = (XX-cogX) * params.shrink_factor + cogX;
  cogY = repmat(cog(2,:),4,1);
  YY = (YY-cogY) * params.shrink_factor + cogY;
  cogZ = repmat(cog(3,:),4,1);
  ZZ = (ZZ-cogZ) * params.shrink_factor + cogZ;
  
  % plot patches
  p = patch(XX,YY,ZZ,params.face_value);
  
% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
