function subseq = subblock_sequence(seq,varargin)
%function subseq = subblock_sequence(seq,ind1[,ind2,...])
%
% For a given cell-array with matrix entries, an identically sized 
% cell-array subseq is generated, in which only the
% matrix-components ind1 x ind2 are copied, e.g. 1:N, 1:M would be
% suitable ranges to be passed in ind1, ind2.

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

subseq = cell(size(seq));
for ci = 1:numel(seq);
  subseq{ci} = seq{ci}(varargin{:});
end;

%| \docupdate 
