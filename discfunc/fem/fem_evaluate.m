function res = fem_evaluate(df, eindices, lcoord, dummy1, dummy2)
%function res = fem_evaluate(df, eindices, lcoord)
%
% method evaluating a fem function df in local coordinates in the point
% `\hat x` = lcoord in several elements simultaneously given by 
% indices eindices.
% res is a length(eindices) x dimrange vector 

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


% dummy arguments is ignored, only required for consistency to 
% other discfunctypes

% Bernard Haasdonk 12.1.2011

% evaluate all reference basis functions in the lcoords, 
% is a ndofs_per_element x dimrange matrix
basis_values = fem_evaluate_basis(df,lcoord);
%res = zeros(params.dimrange,length(eindices));
%res = zeros(length(eindices),params.dimrange);
% linear combination with DOFS to get function values
gid = df.global_dof_index(eindices,:); % only first scalar component!
% insert the incremented dofs
ind = ones(df.dimrange,1)*(1:size(gid,2)); 
ind = ind(:)';
ggid = gid(:,ind);
incr = ones(size(gid,1),1) * (0:(df.dimrange-1)); 
iincr = repmat(incr,1,size(gid,2));
gid = ggid+iincr;
d = df.dofs(gid);
d = reshape(d,length(eindices),df.ndofs_per_element);
res = d * basis_values; 

%fem_evaluate_basis([0,0],params)'*dofs(1,:)'
%| \docupdate 
