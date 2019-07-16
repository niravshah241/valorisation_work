function reduced_data = lin_ds_gen_reduced_data(model,detailed_data)
%function reduced_data = lin_ds_gen_reduced_data(model,detailed_data)
%
% function computing reduced data for rb ds simulation

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


% Bernard Haasdonk 2.4.2009

%model.affine_decomp_mode = 'components';

%W = detailed_data.W;
V = detailed_data.V;
Wtrans = detailed_data.W';
Vtrans = detailed_data.V';

%keyboard;

% reduced system matrices:
%model.t = 0;
x0comp = detailed_data.x0_components;
Q_x0 = length(x0comp);
x0r = cell(1,Q_x0);
for i=1:Q_x0
  x0r{i} = Wtrans * x0comp{i};   
end;

%Acomp = model.A_function_ptr(model);
Acomp = detailed_data.A_components;
Q_A = length(Acomp);
Ar = cell(1,Q_A);
for i=1:Q_A
  Ar{i} = Wtrans * Acomp{i} * detailed_data.V;   
end;

%Bcomp = model.B_function_ptr(model);
Bcomp = detailed_data.B_components;
Q_B = length(Bcomp);
Br = cell(1,Q_B);
for i=1:Q_B
  Br{i} = Wtrans * Bcomp{i};   
end;

%Ccomp = model.C_function_ptr(model);
Ccomp = detailed_data.C_components;
Q_C = length(Ccomp);
Cr = cell(1,Q_C);
for i=1:Q_C
  Cr{i} = Ccomp{i}*detailed_data.V;   
end;

if model.enable_error_estimator
  
  % error estimator matrices:
  M1 = cell(Q_A,Q_A);
  for q1=1:Q_A
  for q2=1:Q_A
    M1{q1,q2} = Vtrans*(Acomp{q1}')*detailed_data.G*Acomp{q2}*detailed_data.V;   
  end;
  end;
  
  M2 = cell(Q_B,Q_B);
  for q1=1:Q_B
    for q2=1:Q_B
      M2{q1,q2} = (Bcomp{q1}')*detailed_data.G*Bcomp{q2};   
    end;
  end;
  
  % M3 is single matrix
  M3 = Vtrans * detailed_data.G * detailed_data.V;
  
  M4 = cell(Q_B,Q_A);
  for q1=1:Q_B
    for q2=1:Q_A
      M4{q1,q2} = (Bcomp{q1}')*detailed_data.G*Acomp{q2}*detailed_data.V;   
    end;
  end;
  
  M5 = cell(1,Q_A);
  for q=1:Q_A
    M5{q} = Vtrans * detailed_data.G * Acomp{q}*detailed_data.V;   
  end;
  
  M6= cell(1,Q_B);
  for q=1:Q_B
    M6{q} = Vtrans * detailed_data.G * Bcomp{q};   
  end;
  
  %m = cell(Q_x0,Q_x0);
  %I_minVW = speye(size(detailed_data.V,1))- detailed_data.V*Wtrans;
  %Mtemp = I_minVW' * detailed_data.G * I_minVW;
  %for q1=1:Q_x0
  %  for q2=1:Q_x0
  %    m{q1,q2} = x0comp{q1}'*Mtemp*x0comp{q2};
  %  end;
  %end;
  
  % instead of large dense matrix allocation equivalent:
  % a' Mtemp b = 
  % a' G b - (a' W) (V' G b) - (a' G V) (W' b) + (a' W) (V' G V) (W' b)
  %
  % these components can not be summed up here,
  % as we must be able to extract online data for other reduced dimensions!!!
  % hence no longer computation of m here 
  % but decomposition as follows:
  %
  % m00: 2d cell array of values x0{q}'G x0{q'} 
  % VtGx0: 1d cell array of of vectors  (V' G x0{q})
  % Wtx0:  1d cell array of of vectors  (W' x0{q})
  % VtGV : matrix V' G V
  %
  % then in rb_simulation the above sum can be easily performed
  
  G = detailed_data.G;
  
  m00 = cell(Q_x0,Q_x0);
  VtGx0 = cell(1,Q_x0);
  Wtx0 = cell(1,Q_x0);
  for q1=1:Q_x0
    v1 = x0comp{q1};
    for q2=1:Q_x0
      v2 = x0comp{q2};
      %    m{q1,q2} = x0comp{q1}'*Mtemp*x0comp{q2};
      %    m{q1,q2} = v1t * G * v2 - (v1t * W) * (Vtrans * G * v2) ...
      %	- (v1t * G * V) * (Wtrans * v2) ...
      %	+ (v1t * W) * (Vtrans * G * V) * (Wtrans * v2);
      m00{q1,q2}= v1'*G*v2;
    end;
    VtGx0{q1} = Vtrans * G * v1;
    Wtx0{q1} = Wtrans * v1;
  end;
  VtGV = Vtrans * G * V;  
end;

reduced_data.x0r = x0r;
reduced_data.Ar = Ar;
reduced_data.Br = Br;
reduced_data.Cr = Cr;
reduced_data.N = size(V,2);

if model.enable_error_estimator
  reduced_data.M1 = M1;
  reduced_data.M2 = M2;
  reduced_data.M3 = M3;
  reduced_data.M4 = M4;
  reduced_data.M5 = M5;
  reduced_data.M6 = M6;
  %reduced_data.m = m;
  reduced_data.m00 = m00; %: 2d cell array of values x0{q}'G x0{q'} 
  reduced_data.VtGx0 = VtGx0; %: 1d cell array of of vectors  (V' G x0{q})
  reduced_data.Wtx0 = Wtx0; %:  1d cell array of of vectors  (W' x0{q})
  reduced_data.VtGV = VtGV; %: matrix V' G V
end;

