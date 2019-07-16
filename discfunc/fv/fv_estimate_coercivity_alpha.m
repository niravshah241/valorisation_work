function [alpha] = fv_estimate_coercivity_alpha(grid,model)
%function [alpha] = fv_estimate_coercivity_alpha(grid,model)
%
% function decomposing L_I = Id + Delta t * M_I and 
% estimating the coercivity constant alpha(mu) of M_I 
% as required for the energy error-estimate in th RB
% simulation.
%
% in case of advdiff should it be the diffusivity constant

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

  
% Bernard Haasdonk 20.7.2006

model.t = 0;
model.dt = model.T / model.nt;

[L_I_diff, bdir_I_diff] = fv_operators_diff_implicit(grid,model);

[L_I_conv, bdir_I_conv] = fv_operators_conv_implicit(grid,model);

[L_I_neu, bneu_I] = fv_operators_neuman_implicit(grid,model);

M_I = L_I_diff + L_I_conv + L_I_neu; % latter two are zero!

% alternative:
%[L_I,L_E,b] = ...
%    operators_conv_diff(grid,model); % independent of U!!     
%M_I = (L_I-speye(size(L_I)))/model.dt;

if model.verbose > 9
  disp('estimating coercivity alpha:')
end;
alpha = abs(eigs(M_I,1,'SM'));
if parmas.verbose > 9
  disp(['alpha = ',num2str(alpha)]);
end;



