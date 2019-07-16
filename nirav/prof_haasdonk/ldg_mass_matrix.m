function M = ldg_mass_matrix(df_info,grid,params)

%function M = ldg_mass_matrix(df_info,grid,params)

%

% function computing the sparse mass matrix of a discrete ldg

% function space on grid with quadrature of degree params.qdeg.

%

%  (M)_ij = \int_Omega phi_i(x) * phi_j(x) dx



% Bernard Haasdonk 31.8.2009

qdeg = params.qdeg;

Mref = ldg_local_mass_matrix(qdeg,df_info);

Mref = reshape(Mref,[size(Mref),1]);

Mseq = repmat(Mref,[1,1,grid.nelements]);

n = size(Mref,1)*grid.nelements;

M = spblkdiag(Mseq); % my own function... 1000 times faster than builtin

diagvals = repmat(grid.A(:)'*2,size(Mref,1),1);

diagvals = diagvals(:);

D = spdiags(diagvals,0,n,n);

M = M * D;



%Mref = sparse(ldg_reference_mass_matrix(df,qdeg));

%n = size(Mref,1)*grid.nelements;

%Mseq = cell(1,grid.nelements);

%[Mseq{:}] = deal(Mref);

%M = blkdiag(Mseq{:});

%diagvals = repmat(grid.A(:)'*2,size(Mref,1),1);

%diagvals = diagvals(:);

%D = spdiags(diagvals,0,n,n);

%M = M * D;%| \docupdate 

