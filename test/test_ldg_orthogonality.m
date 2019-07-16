function OK = test_ldg_orthogonality
%function OK = test_ldg_orthogonality
%
% function testing, whether the ldg basis functions are orthogonal
% or not

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


% Bernard Haasdonk 28.8.2009

OK = 1;

for pdeg = 1:4
  for qdeg = (pdeg*2):11
    [G,V] = ldg_basis_orthonormalization_matrix(pdeg,qdeg);
    maxerr = max(max(abs(G-eye(size(G)))));
    if maxerr>1e-5
      OK = 0;
      disp(['pdeg= ',num2str(pdeg),', qdeg=',num2str(qdeg),...
	    ',|G-I|_infty=',num2str(maxerr)])
      error('ldg functions not orthogonal!')
    end;
  end;
end;%| \docupdate 
