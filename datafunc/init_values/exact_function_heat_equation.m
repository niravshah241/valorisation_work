function res = exact_function_heat_equation(glob, params)

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


t = params.t;

%disp(['time: ', num2str(t)]);

k0 = params.diff_k0;

epsilon = pi;

X = glob(:,1);
Y = glob(:,2);

res = exp(-k0.* epsilon .^2 .* t) ...
       .* sin(1/sqrt(2).*epsilon .* X) ...
       .* cos(1./sqrt(2) .* epsilon .* Y);

