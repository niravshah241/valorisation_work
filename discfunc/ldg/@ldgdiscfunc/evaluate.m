function res = evaluate(df, eindices, lcoord, params)
%function res = evaluate(df, eindices, lcoord, params)
%
% method evaluating a ldg function in local coordinates in the point
% `\hat x = \mbox{lcoord}` in several elements simultaneously given by 
% indices 'eindices'.
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


% Bernard Haasdonk 28.1.2009

if nargin < 4
  params = [];
end
params.ndofs_per_element = df.ndofs_per_element;
params.ndofs = df.ndofs;
res = ldg_evaluate(df.dofs,eindices, lcoord, df.grid, df);

%| \docupdate 
