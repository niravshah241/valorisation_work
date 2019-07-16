function Minv = ldg_inv_mass_matrix(df,grid,qdeg)
%function Minv = ldg_inv_mass_matrix(df,grid,qdeg)
%
% function computing the sparse mass matrix of a discrete ldg
% function on grid with quadrature of degree qdeg.
% df can be a ldgdiscfunc or a structure with fields pdeg and dimrange

% Bernard Haasdonk 31.8.2009

Mref = inv(ldg_local_mass_matrix(qdeg,df));
Mref = reshape(Mref,[size(Mref),1]);
Mseq = repmat(Mref,[1,1,grid.nelements]);

n = size(Mref,1)*grid.nelements;
Minv = spblkdiag(Mseq); % my own function... 1000 times faster than builtin
diagvals = repmat(grid.Ainv(:)'/2,size(Mref,1),1);
diagvals = diagvals(:);
D = spdiags(diagvals,0,n,n);
Minv = Minv * D;%| \docupdate 
