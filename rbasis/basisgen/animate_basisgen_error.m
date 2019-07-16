function animate_basisgen_error(detailed_data, params)

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


save_avi = 1;

params.animation_pause = 0.10;

plot_params.plot_mode = '3d';
plot_params.view      = [-24,22];

switch params.RB_generation_mode
  case 'greedy_uniform_fixed'
    par.numintervals = params.RB_numintervals;
    par.range = params.mu_ranges;
    MMesh = cubegrid(par);
  otherwise
    error('RB_generation_mode not supported for animation');
end

figure;
set(gcf, 'Position', [50 50 800 700]);

if save_avi == 1
  disp('generating basisgen.avi');
  mov = avifile(fullfile(rbmatlabresult,'basisgen_animat.avi'));
  mov.Quality = 100;
  mov = set(mov,'Compression','None');
  frame_rep = ceil(15*params.animation_pause);
end;

errs_sequence = detailed_data.RB_info.post_errs_sequence;
est_sequence  = detailed_data.RB_info.post_est_sequence;

%mu_values     = detailed_data.RB_info.mu_sequence;

vertices =  get(MMesh,'vertex');

for i = 1:size(detailed_data.RB_info.post_errs_sequence,2);

  lambda = 1;%max(1, 10^(floor(log10(...
             %  min(est_sequence(:,i))/(2*max(errs_sequence(:,i))) ...
             % ))));
%  epsilon = 1e-20;
  [err,j] = max(est_sequence(:,i));
  [err,k] = max(errs_sequence(:,i));
%  j = find_vector(mu_values(:,i+1), vertices',epsilon);
  if (isempty(j) || length(j)>1)
    error('vector not or multiply found in vertex list!!');
  end;
  cla;
  plot_leafvertex_data(est_sequence(:,i),MMesh,plot_params); hold on;
  plot_leafvertex_data(errs_sequence(:,i)*lambda,MMesh,plot_params); hold on;
  xlabel('mu_1');
  ylabel('mu_2');
  zlabel(['estimator above error (x', num2str(lambda), ')']);
  zlim([10e-6,10e-1]);
  set(gca,'ZTick',10.^(-5:0));

  pause(params.animation_pause);
  if save_avi
    F = getframe(gcf);
    for fr = 1:frame_rep
      mov = addframe(mov,F);
    end;
  end;
  plot3(vertices(j,1),vertices(j,2),est_sequence(j,i),...
         '.','MarkerSize',30); hold on;
  plot3(vertices(k,1),vertices(k,2),errs_sequence(k,i)*lambda,...
        '.','MarkerSize',30);
  set(gca,'Zscale','log');

  title(['Estimator vs. error plot over M_{train} after RB size N=', ...
         num2str(i), ', at t=T_{end}, with M=', num2str(params.M), ...
         ' and Mstrich=', num2str(params.Mstrich)]);

  if save_avi
    F = getframe(gcf);
    for fr = 1:frame_rep
      mov = addframe(mov,F);
    end;
  end;
  pause(params.animation_pause);
end

if save_avi
  mov=close(mov);
end;

