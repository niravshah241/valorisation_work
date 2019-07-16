function [flux,lambda] = conv_flux_buckley_leverett_derivative(glob, U, params)
% function [flux,lambda] = conv_flux_buckley_leverett_derivative(glob, U, params)
%
% function computing the derivate of a nonlinear convection function for a
% Buckley-Leverett problem.
% ``\partial_{u} f(x,u) = \frac{2M u (1-u)}{(u^2+M(1-U)^2)^2}
% \quad 0\leq u \leq 1``.
%
% required fields of params:
%      bl_k                 : scaling factor for diffusivity `k`
%      bl_M                 : steepness factor of curve `M`
%      U                      : DoF vector of discrete solution
%      grid                   : pointer to grid structure for neighbour
%                               information
%
% See also conv_flux_buckley_leverett()
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

X = glob(:,1);
Y = glob(:,2);

%grid = params.grid;
%
%U = params.U;
%tmpU = repmat(U, 1, grid.nneigh);
%neiU = tmpU;
%real_nb_ind = grid.NBI > 0;
%neiU(real_nb_ind) = U(grid.NBI(real_nb_ind));
%U = 0.5 * (tmpU + neiU);

%vf = @(X,Y) [2*(2*X(:)-pi/2)-(4/3*Y(:)-pi/2)).^2.*cos((2*X(:)-pi/2)...
%             +(4/3*Y(:)-pi/2))+0.1,...
%             2*((2*X(:)-pi/2)-(4/3*Y(:)-pi/2)).^2.*sin((2*X(:)-pi/2)...
%             +(4/3*Y(:)-pi/2))+0.1];
%vf = @(X,Y) [...
%  ((X(:)-Y(:))).^2.*cos((X(:)+Y(:))*pi/2)+0.1, ...
%  ((X(:)-Y(:))).^2.*sin((X(:)+Y(:))*pi/2)+0.1 ...
%  ];

%vf = @(X,Y) [abs((X(:)-Y(:))).^0.50.*cos((X(:)+Y(:))*pi/2)+0.5, ...
%             abs((X(:)-Y(:))).^0.50.*sin((X(:)+Y(:))*pi/2)+0.5];
%vvf = vf(X,Y);
%vvf(X(:)>Y(:),[1,2])=vvf(X(:)>Y(:),[2,1]);


if params.debug && ( max(U)+eps > 1 || min(U) - eps < 0 )
  error('U is outside admissable bounds [0,1]');
end

M = params.bl_M;

%flux = params.bl_k * (2*U(:)*params.bl_M)./(U(:).^2 + params.bl_M).^2;
flux = params.bl_k * ...
       ( 2*M*U(:)    .*   (1-U(:))  ) ...
       ./...
       ( U(:).^2    +  M * (1-U(:)).^2 ).^2;
%flux = [flux, flux];
flux = [ flux, flux ] .* [ones(size(flux)), params.conv_a*X(:)];

lambda = max(flux);

if params.decomp_mode>0
  error('function is nonlinear and does not support affine decomposition!');
end

