function detailed_data=lin_evol_gen_detailed_data(model,model_data)
%function detailed_data=lin_evol_gen_detailed_data(model,model_data)
%
% method, which prepares reduced_data, which is meant as 
% data, that may be dependent on H only, and is not required during the 
% online-simulation, but it may be used for reconstruction purposes or
% computation of online_quantities.  i.e. reduced basis vectors, 
% colateral reduced basis spaces, grid, etc. can be stored here.
%
% So detailed_data is produced and used only in offline-phase algorithms.
%
% allowed dependency of data: H
% allowed dependency of computation: H
% Unknown at this stage: mu, Nmax, N
%
% Required fields of model:
%  gridtype : rectgrid or triagrid to be generated
%  mu_names : cell array of parameter-names
%  mu_ranges : cell array of parameter-intervals
%
% Further fields are required for rb generation,
% see rb_basis_generation
%
% generated fields of detailed_data:
%
%  RB    : reduced basis columns 1,...,Nmax
%  grid  : grid to be used in the subsequent stages
%  RB_info : depending on generation method some detailed information

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


% Bernard Haasdonk 15.5.2007

if(isempty(model_data))
  model_data = gen_model_data(model);
end
detailed_data = [];

skip_rb_generation = false;
detailed_data = structcpy(detailed_data, model_data);
if isfield(model_data, 'mexptr') && isfield(model, 'detailedfname')
  if model_data.mexptr('is_valid_rb')
    if exist(model.detailedfname, 'file')
      tmp = load(model.detailedfname);
      if length(tmp.detailed_data.RB) == model_data.mexptr('get_rb_size')
        disp('loading already generated rb space!');
        detailed_data = tmp.detailed_data;
        skip_rb_generation = true;
      else
        model.mexptr('reset_rb');
      end
    else
      model.mexptr('reset_rb');
    end
  end
end

if ~skip_rb_generation
  detailed_data = rb_basis_generation(model,detailed_data);
  if isfield(model, 'detailedfname')
    save(model.detailedfname, 'model', 'detailed_data');
  end
end


