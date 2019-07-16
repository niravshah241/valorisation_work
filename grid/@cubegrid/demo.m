function demo(dummy)
% function demo(dummy)
% small script demonstrating the possibilities of the cubegrid
% class.

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

% various dimensional grids and plot routine:

disp('');
disp('plotting different dimensional grids and different plot modes')
params.range = {[0,1]};
params.numintervals = [10];
g1 = cubegrid(params);
params.shrink_factor = 0.9;
params.plot_patch = 0;
subplot(1,3,1), plot(g1,params);
title('1D');

params.range = {[0,1],[0,2]};
params.numintervals = [2,4];
g2 = cubegrid(params);
params.plot_patch = 0;
params.axis_equal = 1;
subplot(1,3,2), plot(g2);
title('2D, lines');

params.range = {[0,1],[0,2],[1,2]};
params.numintervals = [2,2,2];
g3 = cubegrid(params);
params.plot_patch = 1;
params.color = [0.5 0.5 0];
params.shrink_factor = 0.9;
subplot(1,3,3), plot(g3, params);
title('3D, patches');

disp('press key to continue');
pause();

% refinement
disp(' ');
disp('demonstration of grid refinement')

% note: indices are leaf-indices, no global element indices!!
g4 = refine(g2,[1,4]);
g4 = refine(g4,[5,6,7]);
g4 = refine(g4,[8,9,10]);

figure;
subplot(1,2,1),plot_leafelement_data(g4,get(g4,'leafelements'),params)
title('element ids of leaf elements');
params.shrink_factor = 1;
subplot(1,2,2),plot_leafvertex_data(g4,1:get(g4,'nvertices'),params) 
title('vertex ids of leaf vertices');

disp('press key to continue');
pause();

disp(' ');
disp('display method of grid:')
display(g4);

disp('press key to continue');
pause();

disp(' ');
disp('get method of grid:')
help cubegrid/get;

end

