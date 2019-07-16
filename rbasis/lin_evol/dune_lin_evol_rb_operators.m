function [LL_I, LL_E, bb, K_II, K_IE, K_EE, m_I, m_E, m ] = ...
      dune_lin_evol_rb_operators(model, detailed_data)

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

  decomp_mode = model.decomp_mode;

  if decomp_mode == 2
    dLL_I = model.coeff_ops.LL_I_ptr(model.mu);
    dLL_E = model.coeff_ops.LL_E_ptr(model.mu);
    dbb   = model.coeff_ops.bb_ptr(model.mu);
%    model.mexptr('set_mu', model.mu);
%    foo   = model.mexptr('rb_operators', decomp_mode);
%    disp('start');
%    dLL_I == foo{1}
%    dLL_E == foo{2}
%    dbb   == foo{3}
%    kron(dLL_I,dLL_E) == foo{6}
%    kron(dLL_E,dbb) == foo{8}
%    kron(dbb,dbb) == foo{9}
    LL_I = [1, dLL_I];
    LL_E = [1, dLL_E];
    bb   = dbb;
    K_II = [1, dLL_I, dLL_I, kron(dLL_I,dLL_I) ];
    K_EE = [1, dLL_E, dLL_E, kron(dLL_E,dLL_E) ];
    K_IE = [1, dLL_I, dLL_E, kron(dLL_I,dLL_E) ];
    m_I  = [dbb, kron(dLL_I,dbb) ];
    m_E  = [dbb, kron(dLL_E,dbb) ];
    m    = kron(dbb,dbb);
  else
    foo   = model.mexptr('rb_operators', decomp_mode);
    dLL_I = foo{1};
    dLL_E = foo{2};
    dbb   = foo{3};
    dK_II = foo{4};
    dK_EE = foo{5};
    dK_IE = foo{6};
    dm_I  = foo{7};
    dm_E  = foo{8};
    dm    = foo{9};

    if decomp_mode == 1

      func2  = @(X) model.dt * X;
      func2t = @(X) model.dt * X';
      func3  = @(X) model.dt.^2 * X;

      identity = speye(size(dLL_I{1}));

      LL_I = [{identity}; ...
        cellfun(func2 , dLL_I, 'UniformOutput', false)];

      LL_E = [{identity}; ...
        cellfun(func2 , dLL_E, 'UniformOutput', false)];

      bb   =  cellfun(func2 , dbb  , 'UniformOutput', false);

      K_II = [{identity}; ...
        cellfun(func2t, dLL_I, 'UniformOutput', false); ...
        cellfun(func2 , dLL_I, 'UniformOutput', false); ...
        cellfun(func3 , dK_II, 'UniformOutput', false)];

      K_EE = [{identity}; ...
        cellfun(func2t, dLL_E, 'UniformOutput', false); ...
        cellfun(func2 , dLL_E, 'UniformOutput', false); ...
        cellfun(func3 , dK_EE, 'UniformOutput', false)];

      K_IE = [{identity}; ...
        cellfun(func2t, dLL_I, 'UniformOutput', false); ...
        cellfun(func2 , dLL_E, 'UniformOutput', false); ...
        cellfun(func3 , dK_IE, 'UniformOutput', false)];

      m_I  = [bb; ...
        cellfun(func3 , dm_I , 'UniformOutput', false)];

      m_E  = [bb; ...
        cellfun(func3 , dm_E , 'UniformOutput', false)];

      m    =  cellfun(func3 , dm   , 'UniformOutput', false);

    else % decomp_mode == 0
      LL_I = 1 + model.dt * dLL_I;
      LL_E = 1 + model.dt * dLL_E;
      bb   =     model.dt * dbb;
      K_II = 1 + model.dt * dLL_I' + model.dt * dLL_I + model.dt^2 * dK_II;
      K_EE = 1 + model.dt * dLL_E' + model.dt * dLL_E + model.dt^2 * dK_EE;
      K_IE = 1 + model.dt * dLL_I' + model.dt * dLL_E + model.dt^2 * dK_IE;
      m_I  =     model.dt * dbb                       + model.dt^2 * dm_I;
      m_E  =     model.dt * dbb                       + model.dt^2 * dm_E;
      m    =                                            model.dt^2 * dm;
    end
  end

%| \docupdate 
