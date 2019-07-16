function reduced_data = lin_stat_gen_reduced_data(model,detailed_data)
%function reduced_data = lin_stat_gen_reduced_data(model,detailed_data)
%
% function computing reduced data.

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

reduced_data = [];

old_mode = model.decomp_mode;

model.decomp_mode = 1; % == components

[A_comp,f_comp] = model.operators(model,detailed_data);

Q_A = length(A_comp);
reduced_data.Q_A = Q_A;
Q_f = length(f_comp);
reduced_data.Q_f = Q_f;

reduced_data.AN_comp = cell(1,length(A_comp));
for q = 1:Q_A
  reduced_data.AN_comp{q}  = ...
      detailed_data.RB'*A_comp{q}*detailed_data.RB;
end;

reduced_data.fN_comp = cell(1,length(f_comp));
for q = 1:Q_f
  reduced_data.fN_comp{q} = ...
      detailed_data.RB' * f_comp{q};
end;

if model.compute_output_functional
  % assumption: nonparametic output functional, then simple RB
  % evaluation possible  
  l_comp = ...
      model.operators_output(model,detailed_data);
  Q_l = length(l_comp);
  reduced_data.lN_comp = cell(1,Q_l);
  for q = 1:Q_l
    reduced_data.lN_comp{q} = detailed_data.RB' * l_comp{q}; 
  end;
end;

N = model.get_rb_size(model,detailed_data);
reduced_data.N = N;

% plus error estimation quantities

% G = (v_r^q, v_r^q) = v_r^q' * K * v_r^q
% with {v_r^q}_q = (v_f^q, v_a^qq')_{qq'}
% G = [Gff, Gfa; Gfa', Gaa];
%
% (v_f^q,v_f^q') = v_f^q' * K * v_f^q

% matrices of coefficient vectors of Riesz-representers:
% K * v_f^q = f^q      (coefficient vector equations)
K = model.get_inner_product_matrix(detailed_data);
% search solution in H10, i.e. set dirichlet DOFs



v_f = zeros(size(detailed_data.RB,1),Q_f);
K_v_f = zeros(size(detailed_data.RB,1),Q_f);
for q = 1:Q_f
  K_v_f(:,q) = f_comp{q};
  v_f(:,q) = K \ f_comp{q};
end;

v_a = cell(N,1);
K_v_a = cell(N,1);
for n = 1:N
  K_v_a{n} = zeros(size(K,1),Q_A);
  v_a{n} = zeros(size(K,1),Q_A);
  for q = 1:Q_A
    K_v_a{n}(:,q) = (A_comp{q}*detailed_data.RB(:,n));
    v_a{n}(:,q) = K \ K_v_a{n}(:,q);
  end;
end;

% compute matrix G components = [G_ff, G_fa{n}; G_fa{n}', G_aa{n}];
reduced_data.Gff = v_f' * K_v_f;
reduced_data.Gfa = cell(1,N);
for n = 1:N
  reduced_data.Gfa{n} = v_f' * K_v_a{n};
end;
reduced_data.Gaa = cell(N,N);
for n1 = 1:N
  for n2 = 1:N
    reduced_data.Gaa{n1,n2} = v_a{n1}' * K_v_a{n2};
  end;
end;

% finally assemble G = [G_ff, G_fa{n}; G_fa{n}', G_aa{n}];
Q_r = Q_f + N * Q_A;
G = zeros(N,N);
G(1:Q_f,1:Q_f) = reduced_data.Gff;
for n = 1:N
  G(1:Q_f,Q_f+(n-1)*Q_A +(1:Q_A)) = ...
      reduced_data.Gfa{n};
  G(Q_f+(n-1)*Q_A +(1:Q_A),1:Q_f) = ...
      reduced_data.Gfa{n}';
end;
for n1 = 1:N
  for n2 = 1:N
    G(Q_f+(n1-1)*Q_A+(1:Q_A),Q_f+(n2-1)*Q_A +(1:Q_A)) = ...
	reduced_data.Gaa{n1,n2};
  end;
end;
reduced_data.G = G;

% set back old model mode
model.decomp_mode = old_mode;
