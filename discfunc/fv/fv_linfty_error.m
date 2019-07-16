function linfty_error = fv_linfty_error(U1,U2,grid,params)
%function linfty_error = fv_linfty_error(U1,U2,[grid],[params])
% compute the infinity-norm error between two Dof vectors.
%
% This function computing the infinity norm error between the two fv-functions
% or function sequences in 'U1,U2': `\|u_1 - u_2\|_{L^\infty(\Omega)}`.
% Actually, the parameters 'grid' and 'params' are currently superfluous, but
% kept for consistency of function line arguments.
%
% parameters:
%  U1 : first vector of Dofs of function `u_1`
%  U2 : second vector of Dofs of function `u_2`
%
% return values:
%  linfty_error: a single value or sequence of values corresponding to the
%                column-differences between 'U1' and 'U2'

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


% Bernard Haasdonk 28.2.2008

%  linfty_error = sqrt(sum(((U1-U2).^2).*repmat(grid.A(:),1,size(U1,2))));
linfty_error = max(abs(U1-U2));

