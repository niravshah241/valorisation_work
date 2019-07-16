function animate_basisgen(detailed_data,params)
%function animate_basisgen(detailed_data,params)
%
% perform an animation of a 2D RB basis construction process for
% grid-based approaches. i.e. the 'RB_generation_mode' in params
% being 'uniform_fixed', 'uniform_refined', 'adaptive_refined'
% the mesh is plotted, the colors indicating the numbers of the
% basis-functions located in the parameter vertices.
% Only 2 parameters, i.e. a 2D parameter grid are plotted currently
%
% optional fields of params:
%      animation_pause : number of seconds after each basis-vector addition 
%                        default is 0.1 sec.

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


% Bernard Haasdonk 29.3.2007

%save_avi = 1; % generate basisgen.avi
save_avi = 0; % generate basisgen.avi

if ~isfield(params,'animation_pause')
  params.animation_pause = 0.05; % pause after each plot
end;

if length(params.mu_names)~=2
  error('currently only 2D parameter spaces are animated!');
end;

params.axis_equal = 0;
mu_values = detailed_data.RB_info.mu_sequence;
mu_frequency = [];
mesh_index_old = 0;
nvertices_old = 0;
vertices_old = zeros(0,2);

switch params.RB_generation_mode
 case {'refined'}
  mesh_index = detailed_data.RB_info.mesh_index_sequence;
  MMesh_sequence = detailed_data.RB_info.MMesh_list;
 case 'uniform_fixed'
  par.numintervals = params.RB_numintervals;
  par.range = params.mu_ranges;
  MMesh0 = cubegrid(par);
  MMesh_sequence = {MMesh0};
  mesh_index = ones(1,size(mu_values,2));
 otherwise
  error('RB_generation_mode not supported for animation.')
end;

figure;

if save_avi == 1
  disp('generating basisgen.avi');
  mov = avifile('basisgen_animat.avi');
  mov.Quality = 100;
  mov = set(mov,'Compression','None');
  frame_rep = ceil(15*params.animation_pause);
end;

for i = 1:length(mesh_index)
  if isempty(find(isnan(mu_values(:,i)),1))
    % get mu_frequency for current collection
    MMesh = MMesh_sequence{mesh_index(i)};
    if (mesh_index(i)> mesh_index_old)
      % extend mu_frequency and check that vertices have remained the same!
      nvertices = get(MMesh_sequence{mesh_index(i)},'nvertices');
      mu_frequency = [mu_frequency, ...
		      zeros(1, nvertices-nvertices_old)];
      vertices =  get(MMesh_sequence{mesh_index(i)},'vertex');
      if ~isequal(vertices(1:nvertices_old,:),vertices_old)
	error('grid vertices changes during refinement!');
      end;
      vertices_old = vertices;
      nvertices_old = nvertices;
    end;
    % find mu_vector in vertex list
    epsilon = 1e-20;
    j = find_vector(mu_values(:,i), vertices',epsilon);
    if (isempty(j) || length(j)>1)
      error('vector not or multiply found in vertex list!!');
    end;
    mu_frequency(j(1)) = mu_frequency(j(1)) + 1;
    cla;
    plot_leafvertex_data(MMesh,mu_frequency, params); hold on;
    p = plot(mu_values(1,i),mu_values(2,i),'.','Markersize',30);
    title([num2str(sum(mu_frequency)),' basis vectors added']);
    pause(params.animation_pause);    
    xlabel(params.mu_names{1});
    ylabel(params.mu_names{2});
   
    if save_avi
      F = getframe(gcf);
      for fr = 1:frame_rep
	mov = addframe(mov,F);
      end;
    end;
  else
    disp(['skipped animation of mu-point ',num2str(i)]);
  end;
end;

if save_avi
  mov=close(mov);
end;

set(p,'visible','off');


% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
