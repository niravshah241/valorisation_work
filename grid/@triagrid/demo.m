function demo(dummy)
% function demo(dummy)
% small script demonstrating the possibilities of the ::triagrid class.
%

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

disp('');
disp('initializing grid by point list and triangles-vertex list');
load('demoinit.mat');
params = [];
g1 = triagrid(p,t,params);

params.shrink_factor = 0.9;
params.plot_patch = 0;
params.color = [1,0,0];
params.axis_equal = 1;
subplot(1,2,1), plot(g1,params);
title('line, shrink');

params.plot_patch = 1;
params.shrink_factor = 1.0;
params.color = [0,1,0];
subplot(1,2,2), plot(g1,params);
title('patch, noshrink');

disp('press key to continue');
pause();

disp('');
disp('elementdata and vertexdata')

figure;
d = sqrt(g1.CX.^2+g1.CY.^2);
subplot(1,2,1);
plot_element_data(g1,d,params);
title('element data');

dv = sin(((g1.X-0.4).^2+(g1.Y-1.0).^2)*10);
subplot(1,2,2);
plot_vertex_data(g1,dv,params);
title('vertex data');

disp('press key to continue');
pause();


disp('');
disp('sequence of data, please move slider')

params.title = 'vertex data sequence';
params.colorbar_location = 'WestOutside'
ndata = 100;
dv = zeros(g1.nvertices,ndata);
for d = 1:ndata
  dv(:,d) = sin(((g1.X-0.4-0.5*d/ndata).^2+(g1.Y-1.0-0.2*d/ndata).^2)*10);
end;
params.plot = @plot_vertex_data;
plot_sequence(dv,g1,params);

disp('press key to continue');
pause();


disp(' ');
disp('demonstration of grid inspect')

inspect(g1);

disp('press key to continue');
pause();

disp(' ');
disp('display method of grid:')
display(g1);

