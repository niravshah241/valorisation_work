function plot_error_landscape(output)
%function plot_error_landscape(output)
% plots an output structure generated by stochastic_error_estimation()
%
% Parameters:
%  output:   output structure generated by stochastic_error_estimation()
%
% Required fields of output:
%  errs:       error values
%  inds:       time index of maximum error values in trajectory
%  run_name:   unique name of error estimation run used for plot titles.
%  bound:      maxmimum error values to be plotted. Above this values error
%              display is cropped.
%  stab_limit: error values below this value are shaded in a different color in
%              surface plots.
%  pf_descr:   plot field description for axis labeling.
%  tsamples:   plot field ranges
%  samples:    combined plot field ranges
%
% Optional fields of output:
%  error_label: label for the y respectively z axis where the error or error
%               estimator is plotted. (default='error')

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


i = logical(output.errs>output.bound);
output.errs(i) = output.bound;
i = logical(isnan(output.errs));
output.errs(i) = output.bound;
C = ones(size(output.errs));
i = logical(output.errs<output.stab_limit);
C(i) = 2;

if ~isfield(output,'error_label')
  output.error_label = 'error';
end

figure;
if length(output.tsamples) == 2
  if isempty(find(C==1,1)); % i.e. all stable
    surf(output.tsamples{1},output.tsamples{2},output.errs');
  else % some not stable
    surf(output.tsamples{1},output.tsamples{2},output.errs',C');
  end;
  shading interp;
  % figure, pcolor(Ms, Ns,C);
  set(gca,'Zscale','log');
  xlabel(output.pf_descr{1});
  ylabel(output.pf_descr{2});
  zlabel(output.error_label);
else
  plot(output.samples, output.errs');
  xlabel(output.pf_descr{1});
  ylabel(output.error_label);
  set(gca,'Yscale','log');
end
run_name = strrep(output.run_name,'_',' ');
title([run_name,' L-infty([0,T],L2) error']);

figure;
if length(output.tsamples) == 2
  surf(output.tsamples{1}, output.tsamples{2}, real(output.inds'));
  xlabel(output.pf_descr{1});
  ylabel(output.pf_descr{2});
  zlabel('time index');
else
  plot(output.samples, output.inds');
  xlabel(output.pf_descr{1});
  zlabel('time index');
end

title('time index of maximum l2-error');
%set(gca,'Zscale','log');
