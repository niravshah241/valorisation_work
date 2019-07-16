function OK = test_power_vector2
%function OK = test_power_vector2
% function testing the powervectors and its derivatives, i.e. check
% whether finite difference approximates the derivative.

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


% Bernard Haasdonk 28.8.2009

OK = 1;

pdeg = 4;

for n = 1:100
  x = rand(1,2);
  h = 1e-7;
  
  for p = 1:pdeg
    D = power_vector2_derivative(x,p);
    
    d1 = (power_vector2(x+[h,0],p)-power_vector2(x,p))/h;
    d2 = (power_vector2(x+[0,h],p)-power_vector2(x,p))/h;
    
    Dappr = [d1,d2]; 
    
    maxerr = max(max(abs(Dappr-D)));
%    disp(['p = ',num2str(p),', maxerr =',num2str(maxerr)]);
    
    if (maxerr>1e-5)
      OK = 0;
      return;
    end;
    
  end;
end;%| \docupdate 
