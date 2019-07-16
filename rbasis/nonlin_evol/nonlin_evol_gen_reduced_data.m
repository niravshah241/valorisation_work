function reduced_data = nonlin_evol_gen_reduced_data(model, detailed_data, params)
%function reduced_data = nonlin_evol_gen_reduced_data(model, detailed_data[, params])
% method which produces reduced_data, which is the data, that will be passed to
% an online-algorithm.
%
% Therefore, no quantities dependent on the high-dimension `H` may be included
% here. Neither may online-data include parameter-dependent `\mu`-quantities.
% So no complete grid or detailed solutions or reduced basis vectors may be
% stored here. So online data is produced in the offline stage, but may be used
% in online-stages. So the computation time may depend on `H`, but the results
% may not depend on this complexity.
%
% allowed dependency of generated data: `N_{\max}, M_{\max}`
% not allowed dependency of data: `H`
% allowed dependency of computation:  `N_{\max}, M_{\max}, H`
% Unknown at this stage: `\mu`,
%
% Generated fields of reduced_data:
%  a0   : a cell array sequence of 'N'-vector components of initial data
%         projection
%  LL_E : a cell array sequence of 'N x N' component-Matrices of explicit
%         operator evaluations.
%  LL_I : a cell array sequence of 'N x N' component-Matrices of implicit
%         operator evaluations.
%  bb_I : a cell array sequence of 'N'-vector components of the implicit offset
%  DE   : 'Nmax x Mmax' cross correlation matrix between 'detailed_data.QM'
%          and 'detailed_data.RB'
%  BM    : 'Mmax x Mmax' interpolation matrix of empirical interpolation.
%  grid_local_ext : part of the grid containing the cells 'TM' plus their
%                   neighbours to be used for local operator evaluation
%  TM_local :       indices of magic-point-elements in the 'grid_local_ext'
%  RB_local_ext :   reduced basis vector values restricted to the support of
%                   'grid_local_ext'
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


% Bernard Haasdonk 16.5.2007

if nargin == 2
  params = [];
end

if ~isfield(model, 'stencil_mode')
  stencil_mode = 'edge';
else
  stencil_mode = model.stencil_mode;
end

if model.verbose > 11
  disp(['stencil mode: ' , model.stencil_mode]);
  disp(['stencil size: ' , num2str(model.local_stencil_size)]);
end

model.decomp_mode = 1;
reduced_data.a0 = model.rb_init_values(model, detailed_data);
Nmax = size(detailed_data.RB,2);
Mmax = size(detailed_data.QM{1},2);

% assuming that components do not change in time, so wlg t = 0!!!
params.t = 0;
[reduced_data.LL_I, reduced_data.bb_I] = ...
    model.rb_operators(model, detailed_data);

A = detailed_data.W;

reduced_data.Nmass = detailed_data.RB' * A * detailed_data.RB;

reduced_data.DE             = cell(size(detailed_data.QM));
reduced_data.BM             = cell(size(detailed_data.QM));
reduced_data.grid_local_ext = cell(size(detailed_data.QM));
reduced_data.RB_local_ext   = cell(size(detailed_data.QM));
reduced_data.TM_local       = cell(size(detailed_data.QM));
reduced_data.Mmass          = cell(size(detailed_data.QM));

for oi = 1:length(detailed_data.QM(:))
  reduced_data.DE{oi} = detailed_data.RB' * A * detailed_data.QM{oi};

  reduced_data.BM{oi} = detailed_data.BM{oi};

  % determine indices of TMmax and its neighbours
  grid = detailed_data.grid;
  TM   = detailed_data.TM{oi};
  if strcmp(stencil_mode, 'vertex') == 1
    mask = zeros(1, grid.nelements);
    nbi  = grid.NBI(TM,:);
    i    = find(nbi > 0);
    % get indices of TM's neighbour-neighbours
    nnbi            = grid.NBI(unique(nbi(i)),:);
    ni              = find(nnbi > 0);
    [nnbuniq,tally] = unique(sort(nnbi(ni)),'first');
    tally           = [tally(2:end);length(nnbi(ni))+1] - tally;
    mask(nnbuniq)   = tally;
    mask(nbi(i))    = 5;
    mask(TM)        = 6;
    mask(mask < 2)  = 0; % skip the elements which have no vertex with the
                         % given elements in TM in common
    eind            = find(mask);
  else
    eind = index_ext(grid,TM,model.local_stencil_size);
  end

  if isa(grid,'rectgrid')|| isa(grid,'triagrid') || isa(grid,'onedgrid')
    reduced_data.grid_local_ext{oi} = gridpart(grid,eind);
%  else % e.g. struct, or not required local grid:
%    disp('please gen nice grid class for 1d grid...')
%    reduced_data.grid_local_ext{oi} = onedgrid_gridpart(grid,eind);  
%    keyboard;
  end;
  reduced_data.RB_local_ext{oi}   = detailed_data.RB(eind,:);
  reduced_data.implicit_crb_index = detailed_data.implicit_crb_index;
  reduced_data.explicit_crb_index = detailed_data.explicit_crb_index;

  % the following messes up the order of the ei-functions!!
  %reduced_data.TM_local = find(mask(eind)==2);

  %goal TM_local(1) is the local element number, which corresponds to
  %the original one, i.e.  eind(TM_local(i)) = TM(i)

  % create kind of 'inversion' map

  % eind : 1:M_local_ext => 1:H

  glob2loc                 = zeros(grid.nelements,1);
  glob2loc(eind)           = 1:length(eind);
  reduced_data.TM_local{oi} = glob2loc(detailed_data.TM{oi});

  % generate grid_local velocity file, if the flux is based on a file:
  % in case of file access, the correct filename must be set here
  % the following is only relevant in case of use of a 
  % params.use_velocitymatrix_file and filecaching mode 2
  if isfield(model,'filecache_velocity_matrixfile_extract') && ...
      (model.filecache_velocity_matrixfile_extract == 2);
    cache_velocity_matrixfile_extract(model, ...
      reduced_data.grid_local_ext{1}.ECX(:,:),...
      reduced_data.grid_local_ext{1}.ECY(:,:),...
      ['M',num2str(model.M),'_','gridpart_Mmax',num2str(Mmax)]);
  end;

  % may be needed by a posteriori error estimator: the CRB mass matrix Q^t W Q
  reduced_data.Mmass{oi} = detailed_data.QM{oi}' * ...
        model.get_inner_product_matrix(detailed_data) * detailed_data.QM{oi};

  % add the localext2global mapping:
  %reduced_data.local_ext_2_global = eind; 

  % check correspondence!
  %disp('Check TM-local-global correspondence. Remove later...');
  %for i = 1:length(detailed_data.TM)
  %  if (eind(reduced_data.TM_local(i)) ~= detailed_data.TM(i))
  %    error('local T_M numbers confused!!!');
  %  end;
  %end;

  % for debugging: add mass-matrix
  %disp('temporary mass matrix is added');
  %reduced_data.M = detailed_data.RB' * A * detailed_data.RB;

end
if(isfield(detailed_data,'time_split_map'))
  reduced_data.time_split_map = detailed_data.time_split_map;
end
reduced_data.N       = model.get_rb_size(model, detailed_data);
reduced_data.M       = cellfun(@(x) size(x,2), detailed_data.QM);
reduced_data.Mstrich = 0;

