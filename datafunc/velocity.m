function vel = velocity(X,Y,params)
%function vel = velocity(X,Y,params)
%
% function computing a vectorial velocity by pointwise
% evaluation in the point sequences indicated by X and Y. 
%
% fields of params:
% name_velocity = 'linear'
%
% in case of 'linear'
% k : scalar diffusivity parameter
%
% Function supports affine decomposition, i.e. different operation modes
% guided by optional field affine_decomp_mode in params. See also the 
% contents.txt for general explanation
%
% optional fields of params:
%   mu_names : names of fields to be regarded as parameters in vector mu
%   affine_decomp_mode: operation mode of the function
%     'none' (default): no parameter dependence or decomposition is 
%                performed. output is as described above.
%     'components': For each output argument a cell array of output
%                 arguments is returned representing the q-th component
%                 independent of the parameters given in mu_names  
%     'coefficients': For each output argument a cell array of output
%                 arguments is returned representing the q-th coefficient
%                 dependent of the parameters given in mu_names  
%
% in 'coefficients' mode, the parameters in brackets are empty

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


% Bernard Haasdonk 11.4.2006

% determine affine_decomposition_mode as integer  
decomp_mode = get_affine_decomp_mode(params);
decomp_mode = params.decomp_mode;
% flag indicating whether the computation respected the decomposition
respected_decomp_mode = 0;

vel = [];

if isequal(params.name_velocity,'linear')
  if decomp_mode == 2
    vel = [params.cx, params.cy];    
  else
    error('no components implemented yet!!');
  end;
  respected_decomp_mode = 1;
else 
  error('velocity function unknown');
end;

if decomp_mode>0 & respected_decomp_mode==0
  error('function does not support affine decomposition!');
end;
 

% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
