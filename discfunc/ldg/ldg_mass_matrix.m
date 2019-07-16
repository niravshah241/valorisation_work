function M = ldg_mass_matrix(df,grid,params)
%function M = ldg_mass_matrix(df,grid,params)
%
% function computing the sparse mass matrix of a discrete ldg
% function on grid with quadrature of degree params.qdeg.
% df can be a ldgdiscfunc or a structure with fields pdeg and dimrange

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


% Bernard Haasdonk 31.8.2009
qdeg = params.qdeg;
Mref = ldg_reference_mass_matrix(df,qdeg);
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
