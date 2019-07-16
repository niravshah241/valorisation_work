function p = plot_sequence(varargin)
%function p = plot_sequence(data,grid,params[,callbackfigure, cbhandle])
%
% plotting a sequence of data slices on polygonal
% 2d grid (constructed from params if empty) and providing a slider
% for changing between the data
% slices. A new figure is opened and the handle returned in p. 
% If further parameters are set, the call is assumed to
% stem from a callback-function
%
% Every column of data is interpreted as one discrete function
% dof vector forwarded to the params.plot() function.
%
% parameters:
%  varargin: usually called with 3 arguments:
%            @code plot_sequence(data, grid, params) @endcode
%   - data - data vector to be plotted
%   - grid - the underlying grid
%   - params - plotting parameters
%   .
%            Alternatively, there can be 5 arguments:
%            @code plot_sequence(data, grid, params, callbackfigure, cbhandle) @endcode
%            This function is set as a callback function for the time slider
%            uicontrol. The additional arguments are:
%   - callbackfigure - handle to the figure
%   - cbhandle       - handle to the object calling the callback function. This
%                      is usually the time slider.
%
% return values:
%  p : figure handle of plot
%
% required fields in params:
%  - 'plot' -- pointer to the plot-function performing the plotting of a
%              single slice, e.g.  plot_element_data(), plot_vertex_data(),
%              fv_plot(), ldg_plot().
%
% optional field of params:
%  - 'title' -- string indicating the title of the newly opened figure
%  - 'clim'  -- 2-vector giving the range of the colormap. If this is set,
%               identical range is used for all slices. Default is the min and
%               max of all slices.
%  - 'clim_tight' -- if this flag is set, the colorbar is set tightly to the
%                    range of every single data slice. Clearly only one of
%                    'clim' or 'clim_tight' should be set.
%
% see also the chosen 'params.plot_function' for its further params-options

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

% open figure
if nargin < 4 % no callback, open new figure
  
%  keyboard;
  f = figure;
  di = 0.1;
  textwidth= 0.15;
  textheight = 0.05;
  txt = uicontrol(f,'Style','text','Units','normalized',...
		  'Position',[di,di,textwidth,textheight],...
		  'String','Slice 1','Tag','text1');
  slider = uicontrol(f,'Style','slider','Units','normalized',...
		     'Position',[di+textwidth+di, di ,1-3*di-textwidth,...
		    textheight],'Tag','slider1',...
		    'Callback','plot_sequence([],[],[],[],gcbf,gcbo)');  
  ax = axes('Units','normalized',...
	    'Position',[di,2*di+textheight,1-2*di,1-3*di - textheight],...
	    'Tag','axes1');
  %  cb = colorbar;
  ud = [];
  ud.data = varargin{1};
  ud.grid = varargin{2};
  ud.params = varargin{3};
  
  if isempty(ud.grid)
    ud.grid = construct_grid(ud.params);
  end;
  
  if isfield(ud.params,'title')
    set(f,'Name',ud.params.title);
  end;

  if ~isfield(ud.params,'clim_tight')
    ud.params.clim_tight = 0;
  end;
  
  if isfield(ud.params,'clim') && ud.params.clim_tight
    error('please only specify one of the parameters clim and clim_tight.');
  end;
  
  if ~isfield(ud.params,'clim') && ~ud.params.clim_tight
    ud.params.clim = [min(min(ud.data)), max(max(ud.data))];
    if ud.params.clim(1) == ud.params.clim(2)
      ud.params.clim = ud.params.clim + [-eps, +eps];
    end;
  end;
  set(f,'Userdata',ud);
  
  sl = findobj(f,'Tag','slider1');
  %set(sl,'min',0);
  %set(sl,'max',ud.params.nt);
  %set(sl,'sliderstep',[1/ud.params.nt,1]);
  set(sl,'min',1);
  if size(ud.data,2)>1
    set(sl,'max',size(ud.data,2));
    set(sl,'sliderstep',[1/(size(ud.data,2)-1),1]);
  else
    set(sl,'max',size(ud.data,2)+eps);
    set(sl,'visible','off');
    %  set(sl,'sliderstep',[0,0]);
  end;
  set(sl,'value',1);
  th = findobj(f,'Tag','text1');
  set(th,'String',('data slice 1') );
  replot(f,sl,1);
  %  set(f,'Resize','on');  
  p = f;
end;

if nargin >=4 % callback-function
  cbf = varargin{4};
  cbo = varargin{5};
  %  disp('performing callback call of plot_fv_data')
  th = findobj(gcbf,'Tag','text1');
  v = round(get(gcbo,'value'));
  set(gcbo,'value',v)
  set(th,'String',['data slice ',num2str(v)]);
  replot(gcbf,gcbo,v);
  p = gcbf;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function replot(fh,oh,v)
ud = get(fh,'Userdata');
%if (size(ud.data,2)<v+1)
%  disp('number of timesteps in U does not match paramt.nt');
%else
%  d = reshape(ud.data(:,v+1),ud.params.nx,ud.params.ny)';

%d = reshape(ud.data(:,v),ud.params.nx,ud.params.ny)';
%xr = ud.params.xrange;
%x = xr(1):(xr(2)-xr(1))/(ud.params.nx):xr(2);
%yr = ud.params.yrange;
%y = yr(1):(yr(2)-yr(1))/(ud.params.ny):yr(2);
%p=pcolor(x,y,[d,zeros(size(d,1),1);zeros(1,size(d,2)+1)]);
%if isfield(ud.params,'no_lines')
%  %    disp('halt in plot_p1_data');
%  %    keyboard;
%  if ud.params.no_lines == 1
%    set(p,'linestyle','None');
%  end;  
%end;
%ax = gca;
%axis equal;
%axis tight;

d = ud.data(:,v);

% if clim_tight, adapt color-range:
if ud.params.clim_tight
  ud.params.clim = [min(min(d)), max(max(d))];
  if ud.params.clim(1) == ud.params.clim(2)
    ud.params.clim = ud.params.clim + [-eps, +eps];
  end;
end;

% delete old data:
ax = findobj(fh,'Tag','axes1');
cla(ax);

%set(ax,'children',[]);

ud.params.plot(ud.grid,d,ud.params);
%feval(ud.params.plot_function,ud.model,ud.grid,d,ud.params);

%disp('after plot...')
%keyboard;

set(gca,'Clim',[ud.params.clim(1), ud.params.clim(2)]);

%ax2 = ax;
%if isfield(ud.params,'show_colorbar') & (ud.params.show_colorbar)
%  if isfield(ud.params,'colorbar_mode') 
%    colorbar(ud.params.colorbar_mode);
%  else
%    colorbar;
%  end;
%  ax2 = gca;
%end;

%if isfield(ud.params,'no_lines') & (ud.params.no_lines)
%  set(p,'LineStyle','None');
%end;
