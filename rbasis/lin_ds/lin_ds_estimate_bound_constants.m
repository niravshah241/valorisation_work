function [C1,C2] = lin_ds_estimate_bound_constants(model,model_data)
%function [C1,C2] = lin_ds_estimate_bound_constants(model,model_data)
%
% function estimating the constant
% ` C1 >= ||expm(A(mu) t)||_G for all t `
%
% and 
%
% ` C2 >= |C(mu)|_G  `

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


% Bernard Haasdonk 3.4.2009

C1 = 0;
C2 = 0;

nmu = model.estimate_lin_ds_nmu;
nX = model.estimate_lin_ds_nX;
nt = model.estimate_lin_ds_nt;

G = model_data.G;

X = rand(size(G,1),nX);
Xnorm = sqrt(sum((G*X).*X,1));
X = X*sparse(1:nX,1:nX,Xnorm.^(-1));
%Xnorm = sqrt(sum((G*X).*X,1));
%disp('check norm 1!')
%keyboard;

M = rand_uniform(nmu,model.mu_ranges);
T = (0:nt)/nt * model.T;

% the following implementation with 2 loops is very bad matlab
% style. But I do not yet see, how that can be vectorized easily...

for mui=1:nmu
  
  tic;
  model = model.set_mu(model,M(:,mui));
  model.decomp_mode = 0; % complete
  model = model.set_time(model,0);
  %model.t = 0; 
  A = model.A_function_ptr(model,model_data);
  C = model.C_function_ptr(model,model_data);

  for ti = 0:nt
%    disp(ti);
    expAtX = expm(A*T(ti+1))*X ;
    normAtX = sqrt(sum((G*expAtX).*expAtX,1));
    
    C1new = max(normAtX);
    if (length(C1new)>1)
      C1new = C1new(1);
    end;
    
    if (C1new>C1) 
      C1 = C1new;
    end;
    
    CX = C*X; 
    if (size(CX,1)>1)
      normCX = sqrt(sum((G*CX).*CX,1));  
    else
      normCX = sqrt(CX.*CX);  
    end;
    
    C2new = max(normCX);
    if (length(C2new)>1)
      C2new = C2new(1);
    end;
    
    if (C2new>C2)
      C2 = C2new;
    end;

    %    if model.debug
    disp(['current C1: ',num2str(C1),' ',...
	  'current C2: ',num2str(C2)]);
    %    end;
    
    t = toc;
    if (ti == 1) && (mui == 1)
      disp(['Estimation of constants C1 and C2 will take ',...
	    num2str(t*(nmu*(nt+1)-2)),' seconds']);
    end;
  end;
    
end;

