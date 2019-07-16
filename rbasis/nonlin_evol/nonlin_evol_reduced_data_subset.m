function reduced_data_subset = nonlin_evol_reduced_data_subset(model,reduced_data)
%function reduced_data_subset = nonlin_evol_reduced_data_subset(model,reduced_data)
% method which modifies reduced_data, which is the data, that will
% be passed to the online-simulation algorithm.
%
% Typically, this routine only does a submatrix extraction of the
% 'reduced_data.M', 'reduced_data.N' sized offline-objects to produce
% 'model.M', 'model.N' sized objects for the real simulation.
%
% Required fields of model:
%  N       : number of reduced basis vectors to choose
%  M       : number of collateral reduced basis vectors to choose for the
%            simulation
%
% Optional fields of model:
%  Mstrich : number of collateral reduced basis vectors to choose for the
%            a posteriori error estimator (default = 0)
%
% Optional fields of reduced_data:
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
%  CE   : N x M interpolation matrix
%  grid_local_ext : part of the grid containing the cells 'TM' plus their
%                   neighbours to be used for local operator evaluation
%  TM_local       : indices of magic-point-elements in the locally extended
%                   grid
%  RB_local_ext   : reduced basis vector values restricted to the support of
%                   'grid_local_ext'
%  s_RB           : in case of a given output functional, this is a vector of
%                   output functional values of the reduced basis
%  s_l2norm       : in case of a given output functional, this is the
%                   `L^2`-norm of the output functional
%  BM    : 'M x M' interpolation matrix of empirical interpolation.
%  DE    : 'N x M' cross correlation matrix between 'detailed_data.QM' and
%          'detailed_data.RB'
%  Mmass : 'M x M' mass matrix computed from 'detailed_data.QM' vectors
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
%  CE   : N x M interpolation matrix
%  grid_local_ext : part of the grid containing the cells TM plus their
%                   neighbours to be used for local operator evaluation
%  TM_local       : indices of magic-point-elements in the locally extended
%                   grid
%  RB_local_ext   : reduced basis vector values restricted to the support of
%                   'grid_local_ext'
%  s_RB           : in case of a given output functional, this is a vector of
%                   output functional values of the reduced basis
%  s_l2norm       : in case of a given output functional, this is the
%                   `L^2`-norm of the output functional
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


% Bernard Haasdonk 16.5.2008

reduced_data_subset = reduced_data;

if ~isfield(model, 'stencil_mode')
  stencil_mode = 'edge';
else
  stencil_mode = model.stencil_mode;
end

if ~isfield(model, 'Mstrich')
  model.Mstrich = 0;
end

if model.N ~= reduced_data.N || any(model.M ~= reduced_data.M) ...
    || model.Mstrich ~= reduced_data.Mstrich

  % extract correct N-sized submatrices and subvectors from reduced_data
  if model.N > length(reduced_data.a0{1})
    error('N too large for current size of reduced basis!');
  end;
  N       = model.N;
  M       = model.M;
  Mstrich = model.Mstrich;
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
  if isfield(reduced_data, 'Nmass')
    reduced_data_subset.Nmass = reduced_data.Nmass(1:N,1:N);
  end

  if isfield(reduced_data, 'BM')

    reduced_data_subset.factor = cell(size(reduced_data.BM));

    for oi = 1:length(reduced_data.BM(:))
      Mmax = length(reduced_data.TM_local{oi});
      MM   = M(oi) + Mstrich;
      if MM > Mmax
        error('M + Mstrich > Mmax');
      end
      reduced_data_subset.DE{oi} = reduced_data.DE{oi}(1:N,1:MM);
      grid = reduced_data.grid_local_ext{oi};

      if strcmp(stencil_mode, 'vertex') == 1
        mask = zeros(1,grid.nelements);
        nbi  = grid.NBI(reduced_data.TM_local{oi}(1:MM),:);
        i    = find(nbi>0);
        mask(nbi(i)) = 2;
        % get indices of TM's neighbour-neighbours
        nnbi            = grid.NBI(unique(nbi(i)),:);
        ni              = find(nnbi > 0);
        [nnbuniq,tally] = unique(sort(nnbi(ni)),'first');
        tally           = [tally(2:end);length(nnbi(ni))+1] - tally;
        mask(nnbuniq)   = tally;
        mask(nbi(i))    = 5;
        mask(reduced_data.TM_local{oi}(1:MM)) = 6;
        mask(mask < 2)  = 0;
%        disp('workaround for Frederikes model')
%        mask(mask < 1)  = 0; % here second neighbours are generated!

        eind            = find(mask);
      else
        % new version with external routine index_ext:
        eind = index_ext(...
          grid,...
          reduced_data.TM_local{oi}(1:MM),...
          model.local_stencil_size);

      end


      reduced_data_subset.grid_local_ext{oi} = gridpart(grid,eind);

      reduced_data_subset.RB_local_ext{oi}   = reduced_data.RB_local_ext{oi}(eind,1:N);

      % the following messes up the order of the ei-functions!!
      %reduced_data.TM_local = find(mask(eind)==2);

      %goal TM_local(1) is the local element number, which corresponds to
      %the original one, i.e.  eind(TM_local(i)) = TM(i)

      % create kind of 'inversion' map

      % eind : 1:M_local_ext => 1:H

      glob2loc       = zeros(grid.nelements,1);
      glob2loc(eind) = 1:length(eind);
      reduced_data_subset.TM_local{oi} = glob2loc(reduced_data.TM_local{oi}(1:MM));
      reduced_data_subset.BM{oi}       = reduced_data.BM{oi}(1:MM,1:MM);
      reduced_data_subset.Mmass{oi}    = reduced_data.Mmass{oi}(1:MM,1:MM);


      % get Q^t W Q for M to M+M'
      if Mstrich > 0
        QWQ  = reduced_data_subset.Mmass{oi}(M(oi)+1,MM);
        reduced_data_subset.factor{oi} = norm(QWQ);
      else
        reduced_data_subset.factor{oi} = NaN;
      end
      % if the velocity file is based on a file, the corresponding
      % subsets must be extracted here:
      % in case of file access, the correct filename must be set here
      % the following is only relevant in case of use of a
      % model.use_velocitymatrix_file and filecaching mode 2
      if isfield(model,'filecache_velocity_matrixfile_extract') && ...
          (model.filecache_velocity_matrixfile_extract == 2);
        % change global velocity file to MMax velocity file generated in
        % offline
        if Mmax ~= model.Mmax
          error('Mmax in model and reduced_data does not correspond!!');
        end;
        model.velocity_matrixfile = ...
          ['gridpart_Mmax',num2str(Mmax),'_',model.velocity_matrixfile];
        % generate M velocity file if not existing
        cache_velocity_matrixfile_extract(model, ...
        reduced_data_subset.grid_local_ext{1}.ECX(:,:),...
        reduced_data_subset.grid_local_ext{1}.ECY(:,:),...
          ['M',num2str(model.M(oi))]);
      end;
    end
  end

  if isfield(model,'name_output_functional')
    reduced_data_subset.s_RB = reduced_data.s_RB(1:N);
    reduced_data_subset.s_l2norm = reduced_data.s_l2norm;
  end;


  reduced_data_subset.N       = N;
  reduced_data_subset.M       = M;
  reduced_data_subset.Mstrich = Mstrich;
  if isfield(reduced_data, 'time_split_map')
    reduced_data_subset.time_split_map = reduced_data.time_split_map;
  end

end

