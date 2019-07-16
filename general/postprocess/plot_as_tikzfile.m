function clim = plot_as_tikzfile(model, params)
% function plot_as_tikzfile(model, params)
% postprocesses a figure and write out an image and a text file that can be
% included in TeX documents.
%
% This method creates three files from a MATLAB figure specified by
% 'params.figure_handle':
%   - A 'png' coded image file containing a snapshot of the figure box
%   - A text file ending '.tikz' that can be included in a LaTeX document
%     adding axes and meta information around the 'png' file.
%   - A text file ending '_test.tex' that can be compiled with pdflatex and
%     outputs a pdf file with the figure output.
%
% Why this is better than including the 'png' file directly:
%  When scaling images for MATLAB figures, often the text for the ticks,
%  legends and axes description is scaled to an unreadable size. Furthermore,
%  text sizes can differ between a LaTeX document and the MATLAB Monospace font
%  type does not look nice in the PDF output. The solution here, circumvents
%  all these problems, because all text used inside the figures will be
%  generated during pdflatex compilation phase.
%
% parameters:
%  params:          Options controlling the output
%
% required fields of params:
%  figure_handle:   handle to the figure that shall be postprocessed
%  filename:        base of filename for the three generated files
%                   'filename.png', 'filename.tikz' and 'filename_test.tex'
%  filepath:        path where the output files shall be stored
%  width:           width of output picture in pixels.
%
% optional fields of params:
%  height:          height of output picture in pixels.
%                   If this parameter does not exist or is set to zero, the
%                   height is calculated from the width and the pictures ratio.
%  save_colorbar:   boolean specifying wether separate files for the colorbar
%                   is generated.
%  print_axes:      boolean specifying wether the tikz file shall include
%                   drawing commands for axes around the picture.
%                   (Default=true)
%  print_axes_label: boolean specifying wether the tikz file shall include
%                   drawing commands for labels at the axes. This field is
%                   ignored if 'print_axes' is set to false. (Default=true)
%  ticks:           integer specifying how many ticks should be added to the
%                   axes drawn around the figure. This field is ignored if
%                   'print_axes' or 'print_axes_label' is set to false.
%                   (Default=5)
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

fn = params.filename;
% dirname
if isfield(params, 'filepath')
  fp = params.filepath;
else
  fp = './';
end
% figure handle
if isfield(params, 'figure_handle')
  fhandle = params.figure_handle;
else
  fhandle = gcf;
end
% colorbar
if isfield(params, 'save_colorbar') && params.save_colorbar
  save_colorbar = 1;
else
  save_colorbar = 0;
end

if ~isfield(params, 'print_axes')
  params.print_axes = 1;
end

if ~isfield(params, 'print_axes_label')
  params.print_axes_label = 1;
end

if ~isfield(params, 'width')
  params.width = 8.0;
end

if ~isfield(params, 'scaleaxis')
  params.scaleaxis = 'scale only axis';
end



if isfield(params, 'ticks')
  ticks = params.ticks;
else
  ticks = 5;
end

% change figure properties for transparencies
set(fhandle,'Color','none');
% get axes:
children = get(fhandle, 'Children');
for i=length(children):-1:1
  if strcmp(get(children(i), 'Type'), 'axes') == 1
    set(children(i),'Color','none');
    ahandle=children(i);
  end
end
%colorbar('peer', ahandle);
%colorbar(ahandle, 'hide');
set(fhandle,'InvertHardCopy','off');
showaxes('hide');

% export figure to png file
im = export_fig(fullfile(fp, fn),'-png','-a3', fhandle);
imsize = size(im);

% save colorbar as well
if save_colorbar
  colorbarfn = fullfile(fp, [fn, 'colorbar']);
  if ~exist([colorbarfn, '.png'], 'file')
    colorh = colorbar('peer',ahandle);
    yticks = get(colorh, 'YTick');
    set(colorh, 'YTick',[]);
    set(colorh, 'YTickLabel',[]);

    cbim = export_fig(colorbarfn, '-png', '-a3', colorh);
%    cbim = export_fig(colorh);
    size(cbim)
  end
end

% color ranges
col_ranges=cell(2,1);
col_ranges{1} = get(gca,'XLim');
col_ranges{2} = get(gca,'YLim');
clim = get(gca,'CLim');

xr = col_ranges{1};
yr = col_ranges{2};

% x-y axes ratio
xyratio = imsize(1)/imsize(2);

% tikz file
fid = fopen(fullfile(fp, [fn, '.tikz']), 'w+');
ftexid = fopen(fullfile(fp, [fn, '_test.tex']), 'w+');

width = params.width;
if ~isfield(params, 'height') || params.height == 0
  height = width * xyratio;
else
  height = params.height;
end

% xtick string
xtstr = '';
ytstr = '';
for i = 1:ticks
  xtstr = [xtstr, ...
           sprintf('%.2fcm/%.1f', (width) * (i-1)/(ticks-1),...
                   (xr(2)-xr(1))*(i-1)/(ticks-1) + xr(1))];
  ytstr = [ytstr, ...
           sprintf('%.2fcm/%.1f', (height) * (i-1)/(ticks-1),...
                   (yr(2)-yr(1))*(i-1)/(ticks-1) + yr(1))];
  if i ~= ticks
    xtstr = [xtstr, ', '];
    ytstr = [ytstr, ', '];
  end
end

if save_colorbar
  colorbarinput = sprintf('   \\\\begin{scope}[xshift=-2cm]\\\\input{%s}\\n\\\\end{scope}',[colorbarfn, '.tikz']);
else
  colorbarinput = '';
end

fprintf(ftexid,['\\documentclass{article}\n',...
    '\\usepackage{tikz,pgfplots}\n',...
    '\\usetikzlibrary{arrows,snakes,shapes,decorations.pathmorphing,backgrounds,fit}\n',...
    '\\begin{document}\n',...
    ' \\begin{tikzpicture}\n',...
    '  \\input{%s}\n',...
    colorbarinput,...
    ' \\end{tikzpicture}\n',...
    '\\end{document}\n'],[fn, '.tikz']);

fclose(ftexid);

if params.print_axes && params.print_axes_label
  ticks = 'data';
else
  ticks = '\empty';
end
fileinput = sprintf(['\\tikzset{background rectangle/.style={fill=white!0}}\n', ...
'\\everymath{\\scriptstyle};\n', ...
'\\begin{axis}[axis on top, enlargelimits=false,\n',...
'xtick=%s,ytick=%s,',params.scaleaxis,', width=%s cm, height=%s cm]\n',...
'  \\addplot graphics[xmin=%s, xmax=%s, ymin=%s, ymax=%s] {%s};\n',...
'\\end{axis}\n'],...
  ticks, ticks, num2str(width), num2str(height), ...
  num2str(xr(1)), num2str(xr(2)), num2str(yr(1)), num2str(yr(2)), fn);

fwrite(fid, fileinput);

fclose(fid);

if save_colorbar
  colorfid = fopen(fullfile(fp, [colorbarfn, '.tikz']), 'w+');
  fileinput = sprintf(['\\tikzset{background rectangle/.style={fill=white!0}}\n', ...
    '\\everymath{\\scriptstyle};\n', ...
    '\\begin{axis}[xtick=\\empty,ytick=%s,enlargelimits=false,\n', ...
    '              height=%s cm,width=50pt]\n', ...
    '  \\addplot graphics[includegraphics={trim=1 1 1 1},xmin=0, xmax=1,\n', ...
    '                     ymin=%s,ymax=%s]\n', ...
    '                    {%s};\n', ...
    '\\end{axis}\n'], ...
    ticks, num2str(height), num2str(clim(1)), num2str(clim(2)), colorbarfn);

  fwrite(colorfid, fileinput);

  fclose(fid);
end

end
