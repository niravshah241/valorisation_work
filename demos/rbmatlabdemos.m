function rbmatlabdemos
%RBMATLABDEMOS  Set up rbmatlab command line demos.
%   This GUI lets you start some common rbmatlab demos.

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


figfile = 'rbmatlabdemos.fig';
if exist(figfile)
  openfig(figfile);
else
  
labelList=str2mat( ...
    'cubical grid demo',...
    'quadratures on edges',...
    'explicit finite volume scheme',...
    'LDG function space',...
    'linear dynamical system',...
    'quadratrues on cells',...
    'RB gui',...
    'Richards evolution problem',...
    'Riemann-Burgers problem',...
    'RB steps',...
    'rectangular grid demo',...
    'Thermalblock problem',...
    'triangular grid demo',...
    'lagrange FE-functions',...
    'FE-interpolation demo',...
    'FEM demo'...
    );

nameList = [...
      'demo_cubegrid          ';
      'demo_edge_quad         ';
      'demo_explicit_FV       ';
      'demo_ldgfunc           ';
      'demo_lin_ds            ';
      'demo_quadratures       ';
      'demo_rb_gui            ';
      'demo_rb_richards_fv    ';
      'demo_rb_riemann_burgers';
      'demo_rb_steps          ';
      'demo_rectgrid          ';
      'demo_thermalblock      ';
      'demo_triagrid          ';
      'demo_femdiscfunc       ';
      'demo_interpol          ';
      'demo_fem               ';
      ];

needswin = zeros(length(namelist),1);
cmdlnwin(labelList, nameList);

end;

