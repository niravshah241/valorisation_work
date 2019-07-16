function reaction = reaction_term(X, Y, U, params)

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


if isfield(params, 'geometry_transformation')

  if isequal(params.geometry_transformation, 'spline')
%    x_hill  = params.geometry_transformation_spline_x;
%    y_hill  = params.geometry_transformation_spline_y;
%    y_hill(2) = params.hill_height;
    p_mu    = spline_select(params);
    [ breaks, coeffs, pieces, order ] = unmkpp(p_mu);
    p_mu_d  = mkpp(breaks, coeffs(:,1:order-1) .* repmat(order-1:-1:1,pieces,1));
    [ breaks, coeffs, pieces, order ] = unmkpp(p_mu_d);
    if order == 1
      p_mu_dd = mkpp(breaks, zeros(size(coeffs)));
    else
      p_mu_dd = mkpp(breaks, coeffs(:,1:order-1) .* repmat(order-1:-1:1,pieces,1));
    end

    div = 1 + ppval(p_mu, X);
    rec = 2 .* (ppval(p_mu_d,X) ./ div).^2 - ppval(p_mu_dd, X) ./ div;

    reaction = 0.001 .* rec .* U;
  end

end;

% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
