function M_array = add_to_array(M,idx,M_array_old)
%function M_array = add_to_array(M,idx,M_array_old)
%
% add M into array M_array_old at position idx
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


if (isempty(M_array_old)) || (length(M_array_old)<idx)
    M_array = M_array_old;
    M_array{idx}=M;
else
    
    M_array = [M_array_old(1:idx-1), M, M_array_old(idx:end)];
%    for i=1:idx-1
%        M_array{i} = M_array_old{i};
%    end
%    M_array{idx}=M;
%    for i=(idx+1):(length(M_array_old)+1)
%        M_array{i}=M_array_old{i-1};
%    end
    
end
    
