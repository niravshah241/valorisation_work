function p = plot_element_data(grid,data,plot_params)
%function p = plot_element_data(grid,data[,plot_params])
% plot method for a 2d polygonal grid and elementwise constant data
% routine can be used for triangular and rectangular grids.
%
% By default a patch plot is performed.
%
% parameters:
%   data : is assumed to be a vector of length 'grid.nelements' with
%          element-wise scalar values
%
% return values:
%   p: is the list of handles to the graphics primitives
%
% optional fields of plot_params:
%   shrink_factor          : if this flag is given, the elements are plotted
%                            shrinked
%   axis_equal             : if this flag is set, set axis to equal scale
%   no_lines               : if this flag is set, no lines are drawn
%   show_colorbar          : if this flag is set, a colorbar is drawn (default
%                            1)
%   colorbar_location      : string specifying the position of the colorbar,
%                            e.g. 'South','EastOutside' (default), etc.
%   clim                   : if this 2-vector is set, the colorbar is set to
%                            these values
%   geometry_transformation: type of transformation function that is to be
%                            applied on the geometry (default 'none')
%   postprocess            : function 'data = f(data,[CX,CY],plot_params)'
%                            postprocessing the given data before it is plotted.
%   colormap               : user-defined colormap
%   plot_type              : string sepcifying wether a 'patch' (default) or a
%                            'contour' plot is created.
%   line_width             : specifyies the line width of plotted lines, e.g.
%                            in a contour plot.
%   contour_lines          : number of contour_lines if plot_type is set to
%                            'contour'. The default is '0' and means that the
%                            number is determined automatically.
%

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


% Bernard Haasdonk 9.5.2007

if nargin<3
  plot_params = [];
end;

if ~(isfield(plot_params,'shrink_factor'))
  plot_params.shrink_factor = 1.0;
end;

if ~(isfield(plot_params,'axis_equal'))
  plot_params.axis_equal = 0;
end;

if ~(isfield(plot_params,'no_lines'))
  plot_params.no_lines = 0;
end;

if ~(isfield(plot_params,'show_colorbar'))
  plot_params.show_colorbar = 1;
end;

if ~(isfield(plot_params,'colorbar_location'))
  plot_params.colorbar_location = 'EastOutside';
end;

if ~(isfield(plot_params,'plot_type'))
  plot_params.plot_type = 'patch';
end

if ~(isfield(plot_params,'contour_lines'))
  plot_params.contour_lines = 0;
end

if ~(isfield(plot_params,'geometry_transformation'))
  plot_params.geometry_transformation = 0;
end;

if (length(data)~=grid.nelements)
  error('length of data does not match number of elements!');
end;

nneigh = grid.nneigh;

% compute vertex coordinates and scale
XX = grid.X(grid.VI(:));
%XX = XX(:);
YY = grid.Y(grid.VI(:));
%YY = YY(:)
CX = grid.CX(:);
CY = grid.CY(:);

DS = grid.DS;

if (plot_params.geometry_transformation)
  [XX,YY] = geometry_transformation(XX,YY,plot_params);
  [CX,CY] = geometry_transformation(CX,CY,plot_params);
end;

if isfield(plot_params, 'postprocess')
  data = plot_params.postprocess(data,[CX,CY],plot_params);
end

if isequal(plot_params.plot_type, 'patch')
  XX = reshape(XX,size(grid.VI)); % nelements*4 matrix
  YY = reshape(YY,size(grid.VI)); % nelements*4 matrix
  CXX = repmat(CX,1,nneigh);
  CYY = repmat(CY,1,nneigh);
  
  % scale coordinates
  XX = (XX - CXX) * plot_params.shrink_factor + CXX;
  YY = (YY - CYY) * plot_params.shrink_factor + CYY;

  % set patch colors for every VERTEX to the value of the function in the CELL
  % center
  CC = repmat(data(:)',nneigh,1);

  p = patch(XX',YY',CC);

  % Change the axes object's background to the background of figure object it
  % the patch is embedded in.
  ax = get(p, 'Parent');
  fi = get(ax, 'Parent');
  set(ax, 'Color', get(fi, 'Color'));

elseif isequal(plot_params.plot_type, 'contour')
  xmin = min(grid.CX);
  xmax = max(grid.CX);
  ymin = min(grid.CY);
  ymax = max(grid.CY);
  lines = sum(grid.CX==xmin);
  cols  = sum(grid.CY==ymin);
  xstep = (xmax-xmin)/(cols-1);
  ystep = (ymax-ymin)/(lines-1);
  if isequal(class(grid), 'rectgrid')
    XX = xmin:xstep:xmax;
    YY = ymin:ystep:ymax;
    cdata = reshape(data,lines,cols)';
  else
    F=TriScatteredInterp(CX,CY,data);
    XX = xmin:xstep/4:xmax;
    YY = ymin:ystep/4:ymax;
    cdata = F(meshgrid(XX,YY));
  end
  if plot_params.contour_lines > 0
    [h,p] = contour(XX,YY,cdata,plot_params.contour_lines);
  else
    [h,p] = contour(XX,YY,cdata);
  end
  line([xmin xmin], [ymin ymin],'Color',[0.9 0.9 0.9]);
  line([xmax xmax], [ymax ymax],'Color',[0.9 0.9 0.9]);
elseif isequal(plot_params.plot_type, 'mesh')
  xmin = min(grid.CX);
  xmax = max(grid.CX);
  ymin = min(grid.CY);
  ymax = max(grid.CY);
  lines = sum(grid.CX==xmin);
  cols  = sum(grid.CY==ymin);
  xstep = (xmax-xmin)/(cols-1);
  ystep = (ymax-ymin)/(lines-1);
  if isequal(class(grid), 'rectgrid')
    XX = xmin:xstep:xmax;
    YY = ymin:ystep:ymax;
    cdata = reshape(data,lines,cols)';
  else
    F=TriScatteredInterp(CX,CY,data);
    XX = xmin:xstep/4:xmax;
    YY = ymin:ystep/4:ymax;
    cdata = F(meshgrid(XX,YY));
  end
  p = surf(XX,YY,cdata);
  set(p, 'edgecolor', 'none');
  cm=colormap([0.7 0.7 0.7; 0.7 0.7 0.7]);
  set(gca,'View',[-2,65]);
  set(gca,'ZLim',[0,1]);
  set(gca,'YLim',[0,1]);
%  cm=colormap('gray');
%  colormap(cm(end:-1:1,:));
  camlight headlight;
  lighting gouraud;

end

if isfield(plot_params,'line_width')
  set(p, 'LineWidth', plot_params.line_width);
end

if plot_params.axis_equal
  axis equal;
  axis tight;
end;

if isequal(plot_params.plot_type, 'patch') && plot_params.no_lines
  set(p,'linestyle','none');
end;

if isfield(plot_params,'clim')
  set(gca,'Clim',plot_params.clim);
end;
if isfield(plot_params,'colormap')
  set(gcf,'Colormap',plot_params.colormap);
end;
if plot_params.show_colorbar
  colorbar(plot_params.colorbar_location);  
end;

