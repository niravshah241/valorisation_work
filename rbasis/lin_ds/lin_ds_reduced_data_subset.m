function reduced_data = lin_ds_reduced_data_subset(model,reduced_data)
%function reduced_data = lin_ds_reduced_data_subset(model,reduced_data)
%
% function extracting reduced data subset for lin_ds reduced
% simulation. The field model.N is the desired new dimension of the
% reduced data and the corresponding submatrices/vectors, etc. of
% reduced_data are extracted.

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

%model.affine_decomp_mode = 'components';

% reduced system matrices:
%model.t = 0;

if model.N== reduced_data.N
  return;  
elseif model.N>reduced_data.N
  error('desired dimensionality N is not covered by online data!!');
else % requested N is smaller => extract subparts 
  N = model.N;
  reduced_data.x0r = subblock_sequence(reduced_data.x0r,1:N);
  reduced_data.Ar = subblock_sequence(reduced_data.Ar,1:N,1:N);
  input_dim = size(reduced_data.Br{1},2);
  reduced_data.Br = subblock_sequence(reduced_data.Br,1:N,1:input_dim);
  output_dim = size(reduced_data.Cr{1},1);
  reduced_data.Cr = subblock_sequence(reduced_data.Cr,1:output_dim,1:N);
  reduced_data.M1 = subblock_sequence(reduced_data.M1,1:N,1:N);
  % reduced_data.M2 = reduced_data.M2; % input_dim does not change
  reduced_data.M3 = reduced_data.M3(1:N,1:N); 
  reduced_data.M4 = subblock_sequence(reduced_data.M4,1:input_dim,1:N);
  reduced_data.M5 = subblock_sequence(reduced_data.M5,1:N,1:N);
  reduced_data.M6 = subblock_sequence(reduced_data.M6,1:N,1:input_dim);
  %reduced_data.m00 = reduced_data.m00; % remains unchanged
  reduced_data.VtGx0 = subblock_sequence(reduced_data.VtGx0,1:N);
  reduced_data.Wtx0 = subblock_sequence(reduced_data.Wtx0,1:N);
  reduced_data.VtGV = reduced_data.VtGV(1:N,1:N);
  reduced_data.N = model.N;
%  keyboard;
end;








