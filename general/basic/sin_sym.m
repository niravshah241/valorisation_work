function Y = sin_sym(X)
%function Y = sin_sym(X)
%
% sine function guaranteeing exact floating-point symmetry
% i.e. sin(x) = sin(x+2*pi) = - sin(-x)
% i.e. usually we get 
%
%  sin(0.1) + sin(0.1+pi) = 4.1633e-17
%  sin_sym(0.1) + sin_sym(0.1+pi) = 4.1633e-17
%  
% Bernard Haasdonk 5.6.2008

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


% averaging over 4 periods:

%generate large lookuptable

resolution = 1000;
step = pi/resolution;
xvals = 0:step:(pi-1e-10);
if length(xvals)~=resolution;
  error('check length of table!!');
end;
sin_tab = sin(xvals);
sin_tab = [sin_tab, 0];
sin_tab(1:end) = 0.5*(sin_tab(1:end)+ sin_tab(end:-1:1));

% map X to interval [0,pi)
X = X(:);
n = floor(X/pi);
odd_n = mod(n,2); % 0 => no flip => 1 flip
X = X - n* pi;
I = round(X/step)+1;
%Y = sin(X);
Y = sin_tab(I);
Y = Y(:);

%keyboard;

%flip sine values 
%Y = Y.*sgns;
fi = find(odd_n);
Y(fi) = -Y(fi);

return;

% old version with averaging, that does not give better results:
%Yplus = sin(X) + sin(X+2*pi)+ sin(X-2*pi)+ X(X+4*pi); 
%Yminus = sin(X-pi) + sin(X-3*pi) + sin(X+pi)+sin(X+3*pi); 
%Y = (Yplus - Yminus)/6; 
%return;

% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
