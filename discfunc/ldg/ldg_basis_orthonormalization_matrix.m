function [G, Vstr] = ldg_basis_orthonormalization_matrix(pdeg,quaddeg)
%function [G, Vstr] = ldg_basis_orthonormalization_matrix(pdeg,quaddeg)
% 
% computation of Gram matrix of ldg basis of degree pdeg with
% quadrature of degree quaddeg. Inner product is L2-inner product
% over reference triangle.
%
% If the resulting matrix is not unity, the strin Vstr defines an
% orthogonal matrix V, that realizes an orthonormalization.
% i.e. if G=X' * X then    (XV)'(XV) is the identity matrix.
% This matrix must be multiplied in evaluate_basis with the
% current basis in order to produce an orthonormal basis.

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


% Bernard Haasdonk 28.1.2009

% set ldg parameters
params.nelements = 1; 
params.pdeg = pdeg;
params.dimrange = 1;

% test orthogonality of basis functions:
G = ldg_local_mass_matrix(quaddeg,params);

%f = @(lcoord) gram_matrix(ldg_evaluate_basis(params,lcoord));
%G = triaquadrature(quaddeg,f);

% determine V such that Id = V' * G * V, i.e. V is the
% coefficient matrix for orthonormalization of the input vectors
G = 0.5*(G+G');
[v,e] = eig(G); % so G = v * e * v',   v' * G * v = e, 
                % e.^(-0.5) * v' * G * v e.^(-0.5) = Id 
V = v * diag(diag(e).^(-0.5));

%keyboard;

if abs(max(V' * G * V -eye(size(G))))>1e-4
  disp('error in orthonomalization matrix')
end;

Vstr = matrix2str(V);
%| \docupdate 
