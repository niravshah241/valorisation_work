function res = lincomb_sequence(seq,sigma)
%function res = lincomb_sequence(seq,sigma)
%
% function performing a linear combination of the elements in the
% cell array seq with coefficients in sigma result is a
% vector/matrix the same size as the entries of seq.
% if sigma = 0, the component is not touched, i.e. the component
% may also be an empty matrix. 

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


% Bernard Haasdonk 15.5.2007

Q = length(sigma);
res = seq{1}* sigma(1);
for q=2:Q
  if sigma(q)~=0
    res = res + sigma(q)*seq{q};
  end;
end;

