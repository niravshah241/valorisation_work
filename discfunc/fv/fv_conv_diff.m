function U = fv_conv_diff(model,model_data)
%function fv_conv_diff(model,model_data)
%
% function performing a full explicit fv-simulation of a
% convection diffusion problem
%
% ``u_t + \nabla \cdot ( f(u,x,t) - k (x,t) \nabla u) = 0``
% ``   u(x,t) = u_0 \quad\mbox{for}\quad t = 0 ``
% ``   u(x,0) = u_{\mbox{dir}} \quad\mbox{for}
%      \quad x \in \partial \Omega_{\mbox{dir}}``
% ``   ( f - k \nabla u) \cdot n = c_{\mbox{neu}}
%      \quad\mbox{for}\quad x \in \partial \Omega_{\mbox{neu}}``
%
% required fields of model:
% T             : final time
% nt            : number of time-intervals until T, i.e. nt+1
%                 solution slices are computed
%
% plus additional parameters required by init_values, conv_flux, and diff_flux

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


% Bernard Haasdonk 25.9.2005

  if nargin ~= 2
    error('wrong number of parameters!');
  end;

  grid = model_data.grid;
  U = zeros(size(grid.CX,1),model.nt+1);
  
  % initial values by midpoint evaluation
  U(:,1) = model.init_values_ptr([grid.CX,grid.CY],model);
  model.dt = model.T/model.nt;
  
  % loop over fv-steps
  for t = 1:model.nt  
    if model.verbose >= 19
      disp(['entered time-loop step ',num2str(t)]);
    end;
    %    if (t==1) || (num_conv_flux.time_dependent)
    if model.verbose >= 19
      disp('computing convective flux ');
    end;
    num_conv_flux = model.num_conv_flux_ptr(model,model_data,U(:,t));
    %    end;
    
    % diffusive flux is always time dependent:
    num_diff_flux = model.num_diff_flux_ptr(model,model_data,U(:,t));

    % Check CFL-condition, Mario Dissertation p. 19: 
    Xi = 0.00001;
    % num_flux.Lg =  max(max(grid.S)) * model.c * 0.25; 
    % Lipschitz-constant of num-flux
    % for LxF = sup_ij |S_ij|* sup_xy |v(x,y)|
    % epsilon = model.k;
    dt_max = (1-Xi)*grid.alpha^3*(grid.hmin^2)/...
             (grid.alpha*num_conv_flux.Lg*grid.hmin+...
              num_diff_flux.epsilon);
    if model.verbose >= 9
      if model.dt < 0.1*dt_max
        if model.verbose > 9
          disp(['currently set dt = ',num2str(model.dt), ...
                ' can be chosen at least 10 times as large.']);
        end;
      elseif model.dt > dt_max
        if model.verbose > 9
          disp(['maximum allowable timestep due to CFL-cond: ',num2str(dt_max)]);
          disp(['currently set dt = ',num2str(model.dt),' is too large!']);
        end;
      end;
    end;

    if model.verbose >= 9
      disp('stop before fv_step')
      if model.verbose >= 100
        keyboard;
      end;
    end;

    NU = fv_conv_diff_step(U(:,t),num_conv_flux,num_diff_flux,grid,model);
    U(:,t+1) = NU;
  end;

  if model.verbose>=5 % linebreak in fprint-list of dots by fv_conv_diff_step
    fprintf('\n');
  end;


