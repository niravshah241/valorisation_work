% script generating a tikz graphic showing trajectories for certain
% selected parameters
%
% required variables that need to be set:
%   - 'model'
%   - 'detailed_data' (including RB information)
%   - 'plot_params'
%   - 'mu_set': a cell array of parameter vectors for which a trajectory
%               shall be generated.
%   - 'imsavepath': a string specifying the directory name, where to put
%                   the generated files to.
%
% optional variables that can be set:
%   - 'colormap'
%   - 'clim'
%   - 'timeinstants': time instants to include into trajectory. Default is
%                     '[ 0, 1/3, 2/3, 1 ] * model.T'
%   - 'filename': file name for the trajectory
%   - 'width': width of single snapshot. Default is 5cm;
%   - 'boxed_snapshots': if set to 'true' snapshots are surrounded by a box

% optimize plot_params
npp                        = plot_params;
if ~isfield(npp, 'plot_type')
  npp.plot_type = 'patch';
end
if isequal(npp.plot_type, 'contour')
  npp.no_lines               = 0;
end
npp.no_axes                = 1;
npp.axes_equal             = 1;
npp.transparent_background = 1;
npp.show_colorbar          = 0;
if exist('colormap','var')
  npp.colormap               = colormap;
end
if exist('clim','var')
  npp.clim                   = clim;
end
if ~exist('timeinstants','var')
  timeinstants = [0, 1/3, 2/3, 1] * model.T;
end
tsteps = max(1,floor(timeinstants/model.T .* model.nt));
if ~exist('filename','var')
  filename = ['trajectory_',model.name];
end
if ~exist('width','var')
  width = 3;
end

% generate reduced_data
reduced_data = gen_reduced_data(model, detailed_data);
model.N = size(detailed_data.RB,2);
model.M = cellfun(@(X)size(X,2), detailed_data.QM,'UniformOutput', true);
model.Mstrich = 0;
model.enable_error_estimator = 0;
reduced_data = extract_reduced_data_subset(model, reduced_data);

subline = find(detailed_data.grid.CY > 0.6 ...
               & detailed_data.grid.CY < 0.6+detailed_data.grid.hmin);
subline_title = {'X'};
subline_values = detailed_data.grid.CX(subline)';

fp = fullfile(imsavepath, params.model_type);
if ~exist(fp, 'dir')
  mkdir(fp);
end

file = '\\pgfplotsset{width=6cm,compat=newest};';
xmulabel = -0.3;
xshift = width+0.2;
yshift = -width-0.2;
ymulabeloff = -0.45*yshift;
ytimelabel = -yshift+0.1;
xtimelabeloff = 0.4*xshift;

for mi=1:length(mu_set)
  newmodel = model.set_mu(model, mu_set{mi});
  newmodel.newton_epsilon=1e-9;
  rb_sim_data = rb_simulation(newmodel, reduced_data);
  rb_sim_data = rb_reconstruction(newmodel, detailed_data, rb_sim_data);

  ymulabel = (mi-1)*yshift+ymulabeloff;
  mu = mu_set{mi};
  muname   = num2str(mu(1));
  for mni=2:length(mu);
    muname = [muname, ', ', num2str(mu(mni))];
  end

  file = [file, sprintf(['\n',...
    '\\draw (%scm,%scm) node [rotate=90]',...
    '{$\\scriptstyle \\boldsymbol \\mu = (%s)$};'], ...
    num2str(xmulabel), num2str(ymulabel), muname)];

  for nt=1:length(tsteps)
    if mi == 1
      xtimelabel = (nt-1)*xshift+xtimelabeloff;
      file = [file, sprintf(['\n',...
        '\\draw (%scm,%scm) node {$\\scriptstyle t = %.3f$};'], ...
        num2str(xtimelabel), num2str(ytimelabel), timeinstants(nt) )];
    end
    npp.plot(model_data.grid, rb_sim_data.U(:,tsteps(nt)), npp);
    set(gca,'PlotBoxAspectRatio',[1 1 1]);

    mustring = strrep(strrep(strrep(strrep(mat2str(mu),' ', '_'),'[',''),']',''),'.','p');
    fn = ['sample_mu_', mustring, '_tstep_', num2str(nt)];

    Uline = rb_sim_data.U(subline,tsteps(nt));
    subline_title = [ subline_title, [mustring,'_t_',num2str(nt) ] ];
    subline_values = [ subline_values; reshape(Uline, 1, length(Uline)) ];

    disp(['Processing ', fn, '...']);
    tikzparams.filename   = fn;
    tikzparams.filepath   = fp;
    tikzparams.width      = width;
    tikzparams.print_axes = 0;
    if(nt == 1 && (mi == 1 || ~exist('clim', 'var')))
      tikzparams.save_colorbar = 1;
      cbfn = [fn,'colorbar'];
    else
      tikzparams.save_colorbar = 0;
    end
    xcoord = (nt-1) * xshift;
    ycoord = (mi-1) * yshift;

    file = [file, sprintf(['\n\n',...
      '\\begin{scope}[xshift=%scm,yshift=%scm]\n',...
      ' \\input{%s}\n'],...
      num2str(xcoord), num2str(ycoord), [fn, '.tikz']),...
      '\end{scope}'];

    tclim = plot_as_tikzfile(model, tikzparams);
    close(gcf);
    pause(1);
  end

end

datafn = 'sublines.dat';
plot(subline_values');
print_datatable(fullfile(fp, datafn), subline_title, subline_values);

% colorbar
ytstr = '';
ticks = 5;
for i = 1:ticks
  ytstr = [ytstr, ...
    sprintf('%.2fcm/%.1f', (-yshift) * (i-1)/(ticks-1),...
    (tclim(2)-tclim(1))*(i-1)/(ticks-1) + tclim(1))];
  if i ~= ticks
    ytstr = [ytstr, ', '];
  end
end
file = [file, sprintf(['\n',...
  '\\begin{scope}[xshift=%scm,yshift=%scm]\n',...
  ' \\pgfdeclareimage[width=%scm,height=%scm]{colorbar}{%s}\n',...
  ' \\pgftext[at=\\pgfpoint{0}{0},left,base]{\\pgfuseimage{colorbar}};\n',...
  ' \\foreach \\y/\\ytext in {%s}\n', ...
  '  \\draw[shift={(0,\\y)}] (1pt,0pt) -- (-1pt,0pt) node[left] (coordsy) {${\\scriptscriptstyle \\ytext}$};\n', ...
  '\\end{scope}\n']',...
  num2str(nt*xshift+1), num2str((mi-1)/2*yshift),...
  num2str(0.2),num2str(-yshift),cbfn,ytstr)];
disp(cbfn);

tfn = fullfile(fp, [filename,'.tikz']);
fid = fopen(tfn, 'w+');
fwrite(fid,file);
fclose(fid);

tfntex = fullfile(fp, [filename,'_out.tex']);
fid = fopen(tfntex, 'w+');
fprintf(fid, ['\\documentclass[a4paper,11pt]{article}\n',...
      '\\usepackage[x11names,rgb]{xcolor}\n',...
      '\\usepackage{tikz,pgfplots}\n',...
      '\\pgfplotsset{plot coordinates/math parser=false}\n',...
      '\\usepackage{amsfonts}\n',...
      '\\usepackage{amsmath}\n',...
      '\\usetikzlibrary{arrows,snakes,shapes,decorations.pathmorphing,backgrounds,fit}\n',...
      '\\usepackage[active,tightpage]{preview}\n',...
      '\\PreviewEnvironment{tikzpicture}\n',...
      '\\oddsidemargin -10mm\n',...
      '\\evensidemargin -10mm\n',...
      '\\topmargin 5mm\n',...
      '\\headheight 0mm\n',...
      '\\headsep 0mm\n',...
      '\\textheight 247mm\n',...
      '\\textwidth 160mm\n',...
      '\\begin{document}\n',...
      ' \\begin{tikzpicture}[scale=0.5]\n',...
      '  \\input{%s}\n',...
      ' \\end{tikzpicture}\n',...
      '\\end{document}\n'],tfn);
fclose(fid);



