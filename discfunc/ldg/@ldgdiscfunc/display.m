function display(df)
% function display(df)
%
% function printing essential information of ldgfunc on screen

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


% Bernard Haasdonk 2.3.2009

disp(['number of elements =',num2str(df.nelements)]);
disp(['dimension of range =',num2str(df.dimrange)]);
disp(['polynomial degree =',num2str(df.pdeg)]);
disp(['number of dofs per element =',num2str(df.ndofs_per_element)]);
disp(['number of dofs =',num2str(df.ndofs)]);
%| \docupdate 
