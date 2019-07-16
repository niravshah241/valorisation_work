function M = rand_uniform(N,intervals)
%function M = rand_uniform(N,intervals)
%
% function generating uniformly distributed random data in a hypercube
% N vectors are generated as columns of the matrix M, intervals is a 
% cell array indicating the borders of the intervals
%  
% example: rand_uniform(100,{[0,1],[100,1000]})

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

  
% Bernard Haasdonk 29.3.2007
  
  mu_dim = length(intervals);
  M = rand(mu_dim,N);
  for i=1:mu_dim
    mu_range = intervals{i};
    M(i,:) = M(i,:)*(mu_range(2)-mu_range(1))+mu_range(1);
  end;  
  
  
% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
