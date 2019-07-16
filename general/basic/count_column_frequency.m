function [unique_vectors, count] = count_column_frequency(vectors)
%function [unique_vectors, count] = count_columns_frequency(vectors)
%
% function removing duplicates from the matrix "vectors" (which
% contains columnwise vectors) and determining the frequency of the
% vectors. i.e. unique_vectors(:,i) appeared count(i) times in the
% original list.
%
% The implementation is "expensive" by a double loop. Could be optimized.

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


% Bernard Haasdonk 26.7.2006

% expensive: double loop!!

nvec = size(vectors,2);

unique_vectors = [];
count = [];

treated = zeros(1,nvec);
for i = 1:nvec;
  if ~treated(i)
    unique_vectors = [unique_vectors, vectors(:,i)];
    % l-infty difference of all other vectors to current one
%    keyboard;
    ii = repmat(vectors(:,i),1,nvec-i+1);
    di = max(abs(vectors(:,i:end)-ii));
    j = find(di==0);
    treated(j+i-1) = 1;
    count = [count, length(j)];
  end;
end;

% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
