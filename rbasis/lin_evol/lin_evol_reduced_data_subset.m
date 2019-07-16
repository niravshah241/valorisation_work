function reduced_data_subset = lin_evol_reduced_data_subset(model,reduced_data)
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
% Optional fields of model:
%  name_output_functional : name of an output functional
%
% Required fields of reduced_data:
%  N    : number of reduced basis vectors in the reduced_data.
%
% optional fields of reduced_data:
%  a0   : a cell array sequence of N-vector components of initial data
%         projection
%  LL_E : a cell array sequence of N x N component-Matrices of explicit
%         operator evaluations.
%  LL_I : a cell array sequence of N x N component-Matrices of implicit
%         operator evaluations.
%  bb   : a cell array sequence of N-vector components of the offset
%  bb_I : a cell array sequence of N-vector components of the implicit offset
%  K_II : a cell array sequence sequence of N x N matrices
%  K_IE : a cell array sequence of N x N matrices
%  K_EE : a cell array sequence of N x N matrices
%  m_I  : a cell array sequence of N-vector components of offset
%  m_E  : a cell array sequence of N-vector components of offset
%  m    : a cell array sequence of scalars (simply copied from 'reduced_data')
%  s_RB     : in case of a given output functional, this is a vector of output
%             functional values of the reduced basis
%  s_l2norm : in case of a given output functional, this is the `L^2`-norm of
%             the output functional
%
% generated fields of reduced_data_subset:
%  a0   : a cell array sequence of N-vector components of initial data
%         projection
%  LL_E : a cell array sequence of N x N component-Matrices of explicit
%         operator evaluations.
%  LL_I : a cell array sequence of N x N component-Matrices of implicit
%         operator evaluations.
%  bb   : a cell array sequence of N-vector components of the offset
%  bb_I : a cell array sequence of N-vector components of the implicit offset
%  K_II : a cell array sequence sequence of N x N matrices
%  K_IE : a cell array sequence of N x N matrices
%  K_EE : a cell array sequence of N x N matrices
%  m_I  : a cell array sequence of N-vector components of offset
%  m_E  : a cell array sequence of N-vector components of offset
%  m    : a cell array sequence of scalars (simply copied from 'reduced_data')
%  s_RB     : in case of a given output functional, this is a vector of output
%             functional values of the reduced basis
%  s_l2norm : in case of a given output functional, this is the `L^2`-norm of
%             the output functional
%  N    : number of reduced basis vectors in the reduced_data.
%
% \note The fields for 'reduced_data_subset' are only generated if they exist in
% 'reduced_data'.

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

reduced_data_subset = reduced_data;

if reduced_data.N ~= model.N
  % extract correct N-sized submatrices and subvectors from reduced_data
  if model.N > length(reduced_data.a0{1})
    error('N too large for current size of reduced basis!');
  end;
  N = model.N;
  if isfield(reduced_data, 'a0')
    reduced_data_subset.a0   = subblock_sequence(reduced_data.a0,1:N);
  end
  if isfield(reduced_data, 'LL_I')
    reduced_data_subset.LL_I = subblock_sequence(reduced_data.LL_I,1:N,1:N);
  end
  if isfield(reduced_data, 'LL_E')
    reduced_data_subset.LL_E = subblock_sequence(reduced_data.LL_E,1:N,1:N);
  end
  if isfield(reduced_data, 'bb')
    reduced_data_subset.bb   = subblock_sequence(reduced_data.bb,1:N);
  end
  if isfield(reduced_data, 'bb_I')
    reduced_data_subset.bb_I = subblock_sequence(reduced_data.bb_I,1:N);
  end
  if isfield(reduced_data, 'K_EE')
    reduced_data_subset.K_EE = subblock_sequence(reduced_data.K_EE,1:N,1:N);
  end
  if isfield(reduced_data, 'K_IE')
    reduced_data_subset.K_IE = subblock_sequence(reduced_data.K_IE,1:N,1:N);
  end
  if isfield(reduced_data, 'K_II')
    reduced_data_subset.K_II = subblock_sequence(reduced_data.K_II,1:N,1:N);
  end
  if isfield(reduced_data, 'm_I')
    reduced_data_subset.m_I  = subblock_sequence(reduced_data.m_I,1:N);
  end
  if isfield(reduced_data, 'm_E')
    reduced_data_subset.m_E  = subblock_sequence(reduced_data.m_E,1:N);
  end
  if isfield(reduced_data, 'm')
    reduced_data_subset.m    = reduced_data.m;
  end
  
  if isfield(reduced_data, 'T')
    reduced_data_subset.T    = reduced_data.T(1:N,1:N);
  end

  if isfield(model,'name_output_functional')
    reduced_data_subset.s_RB = reduced_data.s_RB(1:N);
    reduced_data_subset.s_l2norm = reduced_data.s_l2norm;
  end;

  reduced_data_subset.N = N;

end

