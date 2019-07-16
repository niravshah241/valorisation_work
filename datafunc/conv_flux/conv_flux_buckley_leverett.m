function [flux, lambda] = conv_flux_buckley_leverett(glob, U, params)
% function [flux, lambda] = conv_flux_buckley_leverett(glob, U, params)
% convective flux for Buckley-Leverett problem
%
% function computing the nonlinear convective flux of a Buckley-Leverett
% problem.  ``f(x,u) = \frac{u^2}{u^2+ M (1-u^2)} \quad 0\leq u \leq 1``.
%
% required fields of params:
%      bl_k                    : scaling factor for flux
%      bl_M                    : steepness factor of curve.
%      U                       : DoF vector of discrete solution
%      grid                    : pointer to grid structure for neighbour
%                                information
%
% See also conv_flux_buckley_leverett_derivative().
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


% glob column check
if params.debug
  if ~isempty(glob) && size(glob,1) < size(glob,2)
    warning('coordinates in variable glob are given row-wise, but expected them to be column-wise');
    if params.debug > 2
      keyboard;
    end
  end
end

%grid = params.grid;

%tmpU = repmat(U, 1, grid.nneigh);
%neiU = tmpU;
%real_nb_ind = grid.NBI > 0;
%neiU(real_nb_ind) = U(grid.NBI(real_nb_ind));
%U = 0.5 * (tmpU + neiU);

X = glob(:,1);
Y = glob(:,2);

flux = params.bl_k * U(:).^2./(U(:).^2 + params.bl_M * (1-U(:).^2));
%vf = @(X,Y) [2*((2*X(:)-pi/2)-(4/3*Y(:)-pi/2)).^2.*cos((2*X(:)-pi/2)...
%             +(4/3*Y(:)-pi/2))+0.1,...
%             2*((2*X(:)-pi/2)-(4/3*Y(:)-pi/2)).^2.*sin((2*X(:)-pi/2)...
%             +(4/3*Y(:)-pi/2))+0.1];

%vf = @(X,Y) [...
%  ((X(:)-Y(:))).^2.*cos((X(:)+Y(:))*pi/2)+0.1, ...
%  ((X(:)-Y(:))).^2.*sin((X(:)+Y(:))*pi/2)+0.1 ...
%  ];

%vf = @(X,Y) [abs((X(:)-Y(:))).^0.5.*cos((X(:)+Y(:))*pi/2)+0.5, ...
%             abs((X(:)-Y(:))).^0.5.*sin((X(:)+Y(:))*pi/2)+0.5];
%vvf = vf(X,Y);
%vvf(X(:)>Y(:),[1,2])=vvf(X(:)>Y(:),[2,1]);


flux = [ flux, flux ] .* [ones(size(flux)), params.conv_a*X(:)];

lambda = 1/max(flux(:));

if params.debug && ( max(U)+eps > 1 || min(U) - eps < 0 )
  error('U is outside admissable bounds [0,1]');
end

if params.decomp_mode>0
  error('function is nonlinear and does not support affine decomposition!');
end

