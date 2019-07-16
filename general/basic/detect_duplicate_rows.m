function [un_id, dupl_id, un_id_of_dupl] = ...
    detect_duplicate_rows(matrix, epsilon)
%function [un_id, dupl_id, un_id_of_dupl] = ...
%    detect_duplicate_rows(matrix [, epsilon] );
%
% checks all rows of the matrix on uniqueness by determining, whether
% l2-norm of difference is smaller than epsilon (default 1e-6 if not
% specified) 
% 
% The vector un_id contains all row indices, which belong to unique rows
% and the first occurence of duplicated rows.
% ALl other indices are set into the vector of duplicated row indices 
% dupl_id. 
% The vector un_id_of_dupl gives the indices of the corresponding first
% occurence of the duplicates, e.g.
% the row with index dupl_id(1) is identical to vector un_id_of_dupl(1)
%
% e.g. matrix = [1 2; 2 3; 2 3; 0 0 ; 4 5 ; 2 3] results in
% un_id = [1 2 4 5], dupl_id = [3 6], un_id_of_dupl = [2 2]
%
% A simple loop-version is implemented, which also works for large 
% matrices and rounding errors in large aspect-ratio case are prevented
% by this

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

% Bernard Haasdonk 15.3.2007

if nargin < 2
  epsilon = 1e-6;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% only perform vertorized method up to 1000 rows
%if (size(matrix,1) <= 1000)
%
%% construct l2 distance matrix of all vectors by binomial formula
%
%K1 = repmat(sum(matrix.*matrix,2),1,size(matrix,1));
%K = matrix * matrix';
%Dsqr = K1 + K1' - 2* K;
%dupl_mask = (Dsqr < epsilon^2); 
%
%% a vector is unique, if there is no 1 in its row before the diagonal element
%% So: delete upper triangle matrix including diagonal
%% find true zero rows => uniques.
%% nonzero rows => dupl_ids
%%        take column of first 1 in its row as un_id_of_dupl
%
%[r,c] = ind2sub(size(dupl_mask),find(dupl_mask)); 
%i = find(c>=r);
%rd = r(i);
%cd = c(i);
%s = sub2ind(size(dupl_mask), rd, cd);
%dupl_mask(s) = 0;

%%% test, if upper triangle is really deleted
%%[r,c] = ind2sub(size(dupl_mask),find(dupl_mask)); 
%%i = find(c>=r);
%%if length(i>0)
%%  error('deletion in dupl_mask not successful!');
%%else
%%  disp('deletion in dupl_mask successful, remove test!!');
%%end,

%m = max(dupl_mask,[],2);

%un_id = find(m==0);
%dupl_id = find(m>0);
%[dummy, un_id_of_dupl] = max(dupl_mask(dupl_id,:),[],2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%else % in number of rows larger than 1000, implement trivial loop search
%  warning(['Duplicate rows search by  non-vectorized version for > 1000 ',...
%	   'rows may take some time!!!']);
  
  un_id = [1];
  dupl_id = [];
  un_id_of_dupl  = [];
  
  for r = 2:size(matrix,1)
    % search row r among previous ones
    diff = sum((matrix(1:(r-1),:)-repmat(matrix(r,:),r-1,1)).^2,2);
    diff_mask = (diff < epsilon^2);
    %keyboard;
    [m,i] = max(diff_mask);
    if (m==0) % vector is unique
      un_id = [un_id, r];
    else % vector is duplicate
      dupl_id = [dupl_id, r];
      un_id_of_dupl = [un_id_of_dupl, i];
    end;    
  end;
  
  un_id = un_id(:);
  dupl_id = dupl_id(:);
  un_id_of_dupl = un_id_of_dupl(:);
  
%end;
%keyboard;


% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
