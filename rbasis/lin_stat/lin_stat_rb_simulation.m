function rb_sim_data = lin_stat_rb_simulation(model,reduced_data)
%function rb_sim_data = lin_stat_rb_simulation(model,reduced_data)
% 
% function performing a reduced simulation

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


% B. Haasdonk 22.2.2011

rb_sim_data =[];

old_mode = model.decomp_mode;
model.decomp_mode = 2; % coefficients

[A_coeff,f_coeff] = ...
    model.operators(model,[]);

AN = lincomb_sequence(reduced_data.AN_comp, A_coeff);
fN = lincomb_sequence(reduced_data.fN_comp, f_coeff);

% solution variable:
%uh = femdiscfunc([],model_data.df_info);
%uh.dofs = A\r;    
uN = AN\fN;

% return results:
rb_sim_data.uN = uN;

% plus error estimator
% res_norm = ... % residual norm

% for elliptic compliant case, X-norm (=H10-norm) error estimator:
Q_r = size(reduced_data.G,1);
neg_auN_coeff = -A_coeff * uN'; 
res_coeff = [f_coeff; neg_auN_coeff(:)];
res_norm_sqr = res_coeff' * reduced_data.G * res_coeff;

% direct computation (expensive):
if 0
  % for debugging:
  neg_auN_coeff = neg_auN_coeff(:);
  Q_f = length(f_coeff);
  res_norm_sqr_ff = f_coeff' * reduced_data.G(1:Q_f,1:Q_f) * f_coeff;
  res_norm_sqr_fAu = f_coeff' * reduced_data.G(1:Q_f,(Q_f+1):end) * neg_auN_coeff;
  res_norm_sqr_Auf = neg_auN_coeff' * reduced_data.G((Q_f+1):end,1:Q_f) * f_coeff;
  res_norm_sqr_AuAu = neg_auN_coeff' * reduced_data.G((Q_f+1):end,(Q_f+1):end) * neg_auN_coeff;
  
  model_data = gen_model_data(model);
  detailed_data = gen_detailed_data(model,model_data);
  model.decomp_mode = 0;
  [A,f] = model.operators(model,model_data);
  model.decomp_mode = 2;
  resAu = A * (detailed_data.RB(:,1:model.N) * uN);
  resf = -f;
  res = resAu + resf;
  % residuum functional is res * v 
  % riesz representant: v_r' K v = (v_r,v) = res*v
  % so res = K * v_r;
  K = model.get_inner_product_matrix(detailed_data);
  v_r = K \ res;
  v_rAu = K \ resAu;
  v_rf= K \ resf;
  % res_norm_sqr = (v_r,v_r) = v_r' K v_r = v_r' * res;
  res_norm_sqr2 = v_r' * res;
  res_norm_sqr2_ff = v_rf' * resf;
  res_norm_sqr2_fAu = v_rf' * resAu;
  res_norm_sqr2_Auf = v_rAu' * resf;
  res_norm_sqr2_AuAu = v_rAu' * resAu;
  keyboard;
end;

% prevent possibly negative numerical disturbances:
res_norm_sqr = max(res_norm_sqr,0);    
res_norm = sqrt(res_norm_sqr);
rb_sim_data.Delta = ...
  res_norm/model.coercivity_alpha(model); 
rb_sim_data.Delta_s = ...
  res_norm_sqr/model.coercivity_alpha(model); 

if model.compute_output_functional
  % get operators
  l_coeff = ...
      model.operators_output(model,reduced_data);
  lN = ...
      lincomb_sequence(reduced_data.lN_comp,l_coeff);
  rb_sim_data.s = lN(:)' * rb_sim_data.uN;
  % rb_sim_data.Delta_s = ...;
end;

model.decomp_mode = old_mode;