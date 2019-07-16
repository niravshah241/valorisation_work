function p = plot_discfunc(df,params)
%function p = plot_discfunc(df,params)
%
% function plotting a single scalar discrete function on triangular grid
% A patch plot is performed as default using speified number of
% subsamplings of the triangles. Each patch is plotted linearly 
% interpolated (hence only showing true values in subsampling nodes.)
% On each triangle-edge, subsamp_level points are inserted.
%
% p is the list of handles to the graphics primitives
%
% grid must provide the fields 
%      X,Y,VI,CX,CY,nelements,nvertices and nneigh;
%
% df must specify the discfunc-information: 
%        params.dimrange, params.pdeg
%        evaluate: pointer to local evaluation routine to be called as
%           evaluate(df,1:nel,locs(i,:),grid,params);
%
% optional fields of params:
%   shrink_factor : if this flag is given, the elements are plotted shrinked
%   axis_equal    : if this flag is set, set axis to equal scale
%   no_lines      : if this flag is set, no lines are drawn
%   show_colorbar : if this flag is set, a colorbar is drawn (default 1)
%   colorbar_location : string specifying the position of the
%                 colorbar, e.g. 'South','EastOutside' (default), etc.
%   clim          : if this 2-vector is set, the colorbar is set to
%                   these values
%   subsampling_level : number of intervals per edge

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


% Bernard Haasdonk 2.2.2009

if df.dimrange>1
  error(['plotting of vectorial functions only via extraction of ',...
	' scalar components!'])
end;

if nargin<2
  params = [];
end;

%if ~(isfield(params,'shrink_factor'))
%  params.shrink_factor = 1.0;
%end;

if ~(isfield(params,'axis_equal'))
  params.axis_equal = 0;
end;

if ~(isfield(params,'no_lines'))
  params.no_lines = 1;
end;

if ~(isfield(params,'show_colorbar'))
  params.show_colorbar = 1;
end;

if ~(isfield(params,'colorbar_location'))
  params.colorbar_location = 'EastOutside';
end;

if ~(isfield(params,'subsampling_level'))
  params.subsampling_level = df.pdeg;
end;

%nneigh = grid.nneigh;

% compute vertex coordinates and scale
%XX = grid.X(grid.VI(:));        
%XX = reshape(XX,size(grid.VI)); % nelements*nneigh matrix
%YY = grid.Y(grid.VI(:));       
%YY = reshape(YY,size(grid.VI)); % nelements*nneigh matrix

%CXX = repmat(grid.CX(:),1,nneigh);
%CYY = repmat(grid.CY(:),1,nneigh);

% scale coordinates
%XX = (XX - CXX) *params.shrink_factor + CXX;
%YY = (YY - CYY) *params.shrink_factor + CYY;

%set patch colors 
%CC = data(grid.VI(:));        

% evaluate discrete function

XX_total = zeros(0,size(df.grid.VI,2));
YY_total = zeros(0,size(df.grid.VI,2));
CC_total = zeros(0,size(df.grid.VI,2));

XX = df.grid.X(df.grid.VI(:));        
XX = reshape(XX,size(df.grid.VI)); % nelements*nneigh matrix
YY = df.grid.Y(df.grid.VI(:));       
YY = reshape(YY,size(df.grid.VI)); % nelements*nneigh matrix
nel = df.grid.nelements;
step = 1/(params.subsampling_level+1);

% outer subsampling triangles:
%keyboard;
for i1 = 1:(params.subsampling_level+1);
  for i2 = 1:(params.subsampling_level+2 - i1);
    
    % list of local evaluation points:
    % i-th row is i-th point of the subsamp-triangle
    locs = [ 1-(i1+i2-2) * step,  (i2-1)*step; ...
	     1-(i1+i2-1) * step,   i2 * step; ...
	     1-(i1+i2-1) * step,  (i2-1) * step ];
    % for coordinate computation:
    lincombweights = [locs(:,2), 1 - locs(:,1)-locs(:,2),locs(:,1)];
    % for local evaluating of function
    locs = [1 - locs(:,1)-locs(:,2),locs];
    
    %XX_part = zeros(size(grid.VI));
    %YY_part = zeros(size(grid.VI));
    CC_part = zeros(size(df.grid.VI));
    for i=1:3
      CC_part(:,i) = evaluate(df,1:nel,locs(i,:));
    end;
    XX_part = XX * lincombweights'; 
    YY_part = YY * lincombweights'; 
    
    XX_total = [XX_total; XX_part];
    YY_total = [YY_total; YY_part];
    CC_total = [CC_total; CC_part];
    
  end;
end;

% inner subsampling triangles:
for i1 = 1:(params.subsampling_level);
  for i2 = 1:(params.subsampling_level + 1 - i1);

    % list of local evaluation points:
    % i-th column is i-th point of the subsamp-triangle
    locs = [ 1-(i1+i2-1) * step,(i2-1)*step; ...
	     1-(i1+i2-1) * step, i2   * step; ...
	     1-(i1+i2) * step  , i2 * step ];
    lincombweights = [locs(:,2), 1 - locs(:,1)-locs(:,2),locs(:,1)];
    locs = [1 - locs(:,1)-locs(:,2),locs];

    %XX_part = zeros(size(grid.VI));
    %YY_part = zeros(size(grid.VI));
    CC_part = zeros(size(df.grid.VI));
    for i=1:3
      CC_part(:,i) = evaluate(df,1:nel,locs(i,:));
    end;
    XX_part = XX * lincombweights'; 
    YY_part = YY * lincombweights'; 

    XX_total = [XX_total; XX_part];
    YY_total = [YY_total; YY_part];
    CC_total = [CC_total; CC_part];
    
  end;
end;

%different options with different results ???
%figure,
p = patch(XX_total',YY_total',CC_total');
%figure,
%p = patch(XX_total(1:10,:)',YY_total(1:10,:)',CC_total(1:10,:)');
%hold on;
%%p = patch(XX_total(11:end,:)',YY_total(11:end,:)',CC_total(11:end,:)');
%keyboard;

if params.axis_equal
  axis equal;
  axis tight;
end;

if params.no_lines
  set(p,'linestyle','none');
end;

%keyboard;

if params.show_colorbar
  if isfield(params,'clim')
    set(gca,'Clim',params.clim)
  end;
  colorbar(params.colorbar_location);  
end;
hold on;

if ~params.no_lines
  p2 = plot(df.grid);
  set(p2,'Color','k')
  p = [p(:);p2];
end;

%| \docupdate 
