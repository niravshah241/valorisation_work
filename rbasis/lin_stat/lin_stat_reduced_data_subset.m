function reduced_data_subset = lin_stat_reduced_data_subset(model,reduced_data)
%function reduced_data_subset = lin_evol_reduced_data_subset(model,reduced_data)
% method which modifies reduced_data, which is the data, that will
% be passed to the online-simulation algorithm.
%
% Typically, this routine only does a submatrix extraction of the 'Nmax' sized
% offline-objects to produce 'N' sized objects for the real simulation.
%
% Required fields of model:
%  N       : number of reduced basis vectors to choose
%
% Required fields of reduced_data:
%  N    : number of reduced basis vectors in the reduced_data,
%         must be larger than model.N!!

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


% Bernard Haasdonk 27.2.2011

reduced_data_subset = reduced_data;

if reduced_data.N ~= model.N
  % extract correct N-sized submatrices and subvectors from reduced_data
  if model.N > reduced_data.N
    error('N too large for current size of reduced basis!');
  end;
  N = model.N;
  if isfield(reduced_data, 'AN_comp')
    reduced_data_subset.AN_comp   = ...
	subblock_sequence(reduced_data.AN_comp,1:N,1:N);
  end
  if isfield(reduced_data, 'fN_comp')
    reduced_data_subset.fN_comp = subblock_sequence(reduced_data.fN_comp,1:N);
  end
  if isfield(reduced_data, 'lN_comp')
    reduced_data_subset.lN_comp = subblock_sequence(reduced_data.lN_comp,1:N);
  end
  % inner product matrix of riesz-representers of residual components:
  reduced_data_subset.Gff = reduced_data.Gff;
  Q_f = length(reduced_data_subset.fN_comp);
  Q_a = length(reduced_data_subset.AN_comp);
  reduced_data_subset.G = reduced_data.G(1:(Q_f+N*Q_a),1:(Q_f+N*Q_a));
  
  reduced_data_subset.N = N;

end

