function B = model_PCA_fixspace(model,model_data,X,XFix,k)
%function B = model_PCA_fixspace(model,model_data,X,XFix [, [A], [k] ])
% 
% function computing a PCA basis of the set of column vectors X projected
% on the orthogonal complement of span(XFix)
% The computation is performed explicitly using the 
% gram matrix.
% Additionally, the sparse matrix of the inner product A on the
% vector space can be provided, i.e.
% i.e. <x,x'> =  x'*A*x
% Additionally, the number of principal components can be specified

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


% 5.10.2005 Bernard Haasdonk

%epsilon = 1e-10;
% epsilon = eps; % = 2e-16

% trivial weighting if matrix is not given:

A = model_data.W;

if nargin < 4
  k = size(X,2);
end;

% the following is expensive, use of matlab-function instead ?
XFixON = model_orthonormalize(model,model_data,XFix);
%XFixON = delzerocolumns(XFixON,[],A);

if ~isempty(XFix)
     Xo = X - XFixON * (XFixON' * A * X); 
else 
    Xo = X;
end;

%if size(Xo,1)<size(Xo,2) % ordinary correlation matrix
%  [e,v] = eig(Xo*Xo');
%  e = e(:,end:-1:1);
%  evalues = diag(v);
%  evalues = evalues(end:-1:1);
%  fi = find(abs(evalues)>epsilon);
%  e = e(:,fi);
%else  
% via gram matrix of orthogonalized trajectory
K = Xo'*A*Xo;
K = 0.5*(K+K'); % required for rounding problems

% generate descending list of eigenvectors/values:
if k<size(K,2)
  [ep,vp] = eigs(K,k);
  evalues = diag(vp);
else
  [ep,vp] = eig(K);
  ep = ep(:,end:-1:1);
  evalues = diag(vp);
  evalues = evalues(end:-1:1);
end;
%fi = find(abs(evalues)>=epsilon);
fi = 1:length(evalues);

% project Xo vectors on eigenvectors
B = Xo * ep(:,fi) * diag(evalues(fi).^(-0.5));
%end;

% ensure that only real valued vectors are returned
while (~isreal(B))
%  disp('complex eigenvector occured: please check!');
%  keyboard;
  B= B(:,1:end-1);
end;

% the following is theoretically superfluous, but does improve the 
% orthogonality, i.e. e' A e => identity up to 1e-8
%                     B' A B => identity up to 1e-16
%B = orthonormalize(B,A);

%| \docupdate 
