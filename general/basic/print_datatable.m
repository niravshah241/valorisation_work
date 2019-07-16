function print_datatable(fname, title, values, blocksize)
% function print_datatable(fname, title, values, blocksize)
%
% print a matrix in a tabulator separated table stored in file 'fname'. This
% can be used by all pgfplots commandos in LaTeX.
%
% arguments:
% 'fname': file name of file where the table is stored
% 'title': vector of column titles
% 'values': matrix whose rows are written to the table
% 'blocksize': separate blocksize values with a new line in each column,
% needed for patch plots
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


if nargin < 4
  blocksize = size(values,2);
end

numblocks = size(values,2)/blocksize;

assert(numblocks-floor(numblocks)==0);

fid=fopen(fname,'w');

assert(length(title) == size(values, 1));

fpstring  = [repmat('%s\t', 1, length(title)), '\n'];
fprintf(fid, fpstring, title{:});
fpstring2 = [repmat('%e\t', 1, length(title)), '\n'];
for i=1:numblocks
  fprintf(fid, fpstring2, values(:,(i-1)*blocksize+(1:blocksize)));
  fprintf(fid, '\n');
end

fclose(fid);
