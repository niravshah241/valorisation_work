function reduced_data = rb_output_functional_prep(...
    model,reduced_data,detailed_data)
%function reduced_data = rb_output_functional_prep(...
%    model,reduced_data,detailed_data)
%
% function computing the output-functional evaluation of the
% reduced basis RB in detailed_data into the aditional fields s_RB
% and s_l2norm.
%
% Required fields of model:
% name_output_functional  : name of output functional to compute

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

% Bernard Haasdonk 16.5.2008

nRB = size(detailed_data.RB,2);
s_RB = zeros(nRB,1);

for n = 1:size(detailed_data.RB,2)
  s = output_functional(model,detailed_data,detailed_data.RB(:,n));
  s_RB(n) = s.value;
end;
reduced_data.s_RB = s_RB;
reduced_data.s_l2norm = s.l2norm;

