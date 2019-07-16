function h=nonlin_evol_plot_detailed_data(model, detailed_data, plot_params)
%function h=nonlin_evol_plot_detailed_data(model, detailed_data, plot_params)
% plot the reduced basis, the colateral reduced basis, the interpolation
% points and the maximum error decrease during CRB generation.
%
% return values:
%   h : figure handle of plot

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


% Bernard Haasdonk 23.5.2007

%plot_RB_sequence(detailed_data,model);



if ~isfield(plot_params,'plot')
  plot_params.plot = model.plot;
end;

if isfield(detailed_data,'RB') && ~isempty(detailed_data.RB)
  tpp = plot_params;
  tpp.clim = [];
  tpp.title = 'Orthonormal reduced basis';
  tpp = rmfield(tpp, 'clim');
  tpp.clim_tight = 1;
  plot_sequence(detailed_data.RB, detailed_data.grid, tpp);
else
  disp('warning: no RB in detailed data!! not suitable for simulation yet!');
end;

if isfield(detailed_data, 'QM')

num_interpolations = min(3,length(detailed_data.QM));

for i = 1:num_interpolations
  tpp = plot_params;
  tpp.title = ['Interpolation basis functions q_m for operator no. ', ...
               num2str(i)];
           tpp.clim = [];
  tpp = rmfield(tpp, 'clim');
  tpp.clim_tight = 1;
  plot_sequence(detailed_data.QM{i}, detailed_data.grid, tpp);
  figure;
end

disp(['Interpolation matrices B should have lower-triangular structure,' ...
      'one-diagonal:']);
for i = 1:num_interpolations
  BM_n = length(detailed_data.BM{i});
  mask_upper_triangular_block = repmat(1:BM_n,BM_n,1)>repmat((1:BM_n)',1,BM_n);
  if(all(abs(detailed_data.BM{i}(mask_upper_triangular_block)) < 10000*eps))
    disp([num2str(i), 'th matrix B fulfills the condition']);
  else
    if model.debug
      disp('The following matrix is not of lower-triangular form');
      pause
      disp(detailed_data.BM{i});
    end
  end
end

for i = 1:num_interpolations
  u = zeros(detailed_data.grid.nelements,1);
  u(detailed_data.TM{i}) = 1;
  plot_params.plot(detailed_data.grid, u, plot_params);
  title(['Interpolation points/DOFS for operator no.', num2str(i)]);

  figure;
end

% plot selected snaphots (mu_1, mu_2, t^k) in a 3d plot in case of 2d parameter
% space
if length(model.mu_names)==2
  for i = 1:num_interpolations
    Msize  = size(detailed_data.QM{i},2);
    colorm = colormap(hot(Msize));
    colorm = colorm(Msize:-1:1,:);
    coord = [detailed_data.ei_info{i}.extension_mus(:,1:Msize); ...
             detailed_data.ei_info{i}.extension_filepos(1:Msize)'];
    scatter3(coord(1,:),coord(2,:),coord(3,:),20,colorm,'filled');
    xlabel(['\mu_1 = ',model.mu_names{1}]);
    ylabel(['\mu_2 = ',model.mu_names{2}]);
    %ylabel(['\mu_2 = c_{init}']);%,model.mu_names{2}]);
    zlabel('time index k');
    title(['snapshots selected for collateral basis for operator no.', num2str(i)])
    set(gca,'Box','on');
    figure;
  end
else
  disp('skipped plot of ei-snapshot-distribution, dimension of mu not 2');
end;

h=[];
for i = 1:num_interpolations
  h=plot(detailed_data.ei_info{i}.max_err_sequence);
  set(gca,'Yscale','log');
  title(['EI-interpol error decrease for operator no.', num2str(i)]);
  if i < num_interpolations
    figure;
  end
end

end
