function display(grid)
%function display(grid)
% display method for ::triagrid
%
% This inherits ::gridbase.display() and adds information on
%   - the number of interior edges and
%   - the number of boundary edges.

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


% Bernard Haasdonk 9.5.2007

% simple forward to gridbase/display
display@gridbase(grid);
disp(['number of interior edges  : ', num2str(grid.nedges_interior)]);
disp(['number of boundary edges  : ', num2str(grid.nedges_boundary)]);

