function detailed_data = lin_ds_gen_detailed_data(model,model_data)
%function detailed_data = lin_ds_gen_detailed_data(model,model_data)
%
% function computing the detailed data for a ds reduced basis simulation
% i.e. the reduced basis, matrix G, components of A,B,C, etc.
%
% Generated fields of detailed_data:
% V:   reduced basis
% W:   reduced projection basis (biorthonormal to V)

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


% Bernard Haasdonk 2.4.2009

if ~model.affinely_decomposed
  error(['rb simulation only possible for affine parameter dependent ',...
	 'linear dynamical systems']);
end;

%detailed_data = [];
% field RB is generated
detailed_data = rb_basis_generation(model,model_data);

%if isfield(model,'RB_filename')
%  tmp = load(fullfile(rbmatlabresult,model.RB_filename));
%  detailed_data.V = tmp.V;
%  detailed_data.W = tmp.W;
%else
%  detailed_data.V = eye(3,2); 
%  detailed_data.W = eye(3,2); % biorthonormal 
%end;
%detailed_data.RB_info = tmp.RB_info;
%detailed_data.V = tmp.RB;
%detailed_data.W = detailed_data.V; % biorthogonal as V orthogonal
detailed_data.G = model_data.G;
detailed_data.x0_components = model_data.x0_components;
detailed_data.A_components = model_data.A_components;
detailed_data.B_components = model_data.B_components;
detailed_data.C_components = model_data.C_components;



























