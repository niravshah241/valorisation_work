function success = test_affine_decomp(fct,nargout,params_pos,varargin)
%function success = test_affine_decomp(fct,nargout,params_pos,varargin)
%
% function calling the function given by handle fct with the
% arguments varargin and expecting nargout many output arguments.
% The fct is called with the 'varargin{params_pos}.decomp_mode' set to 1,2 and 
% linear combination is performed and checked, if result is
% identical to the 'complete' mode 'decomp_mode =0'.
%
% e.g. the routine
%
% @code
%      [A1, A2, A3] = my_fct(params,data)
% @endcode
%
% is called with
%
% @code
%    test_affine_decomp(@my_fct,3,1,params,data);
% @endcode
%
% where params.decomp_mode is expected to be 0.

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


% Bernard Haasdonk 7.9.2009

if varargin{params_pos}.decomp_mode ~=0
  error('please call this routine only in decomp_mode == 0 modus!');
end;

var_complete = cell(1,nargout);
%[var_complete{1:nargout}]= fct(params,varargin{:});
[var_complete{1:nargout}]= fct(varargin{:});

varargin{params_pos}.decomp_mode = 2; 
%params.decomp_mode = 2;
var_coefficients = cell(1,nargout);
[var_coefficients{1:nargout}]= fct(varargin{:});

%params.decomp_mode = 1;
varargin{params_pos}.decomp_mode = 1; 
var_components = cell(1,nargout);
[var_components{1:nargout}]= fct(varargin{:});

var_assembled = cell(1,nargout);
err = zeros(1,nargout);
success = 0;
for q = 1:nargout
  if isempty(var_components{q})
    if ~isempty(var_coefficients{q})
      error('components empty but not coefficients!!');
    end;
  else
    var_assembled{q} = lincomb_sequence(var_components{q},var_coefficients{q});
    err(q) = max(abs(var_assembled{q}(:) - var_complete{q}(:)));
    if (err(q)>1e-8)
      error(['error in affine decomposition of return argument ',...
	     num2str(q),'!']);
    end;
  end;
end;

success = 1;
%| \docupdate 
