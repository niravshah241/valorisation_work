function res = lincomb_sequence3(seq,sigma1,sigma2,sigma3)
%function res = lincomb_sequence3(seq,sigma1,sigma2,sigma3)
%
% function performing a linear combination of the elements in the
% 3d cell array seq with coefficients in sigma1, sigma2, sigma3 result is a
% vector/matrix the same size as the entries of seq.
% size of seq is length(sigma1) x length(sigma2) x length(sigma3)
% if some sigma = 0, the component is not touched, i.e. the component
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

  Q1 = length(sigma1);
  Q2 = length(sigma2);
  Q3 = length(sigma3);
  res = zeros(size(seq{1,1,1}));
  for q1=1:Q1
    if sigma1(q1)~=0 
      for q2=1:Q2
	if sigma2(q2)~=0 
	  for q3=1:Q3
	    if sigma3(q3)~=0      
	      res = res + sigma1(q1)*sigma2(q2)*sigma3(q3)* seq{q1,q2,q3};
	    end;
	  end;
	end;
      end; 
    end;
  end;
