function V = fem_basis_weight_matrix(pdeg)
%function V = fem_basis_weight_matrix(pdeg)
%
% function returning the weight matrix for a lagrange basis on the
% reference triangle used by fem_evaluate_basis
% l(x) = V * p(x)
% with p(x) the monomial basis, i.e. power_vector2
%
% for pdeg = 1 the matrix is explicitly given
% for higher pdegs it is computed, hence expensive. Simply insert
% the obtained matrices hardcode into this file, if your pdeg is
% not yet available.

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


% Bernard Haasdonk 12.1.2011

switch pdeg
 case 1
   V = [1,-1,-1 ; ...
	0, 1, 0; ...
	0, 0, 1];  
 case 2
  V = [1    -3    -3     2     4     2
     0     4     0    -4    -4     0
     0    -1     0     2     0     0
     0     0     4     0    -4    -4
     0     0     0     0     4     0
     0     0    -1     0     0     2
      ];
 case 3
  V = [
     2   -11   -11    18    36    18    -9   -27   -27    -9
     0    18     0   -45   -45     0    27    54    27     0
     0    -9     0    36     9     0   -27   -27     0     0
     0     2     0    -9     0     0     9     0     0     0
     0     0    18     0   -45   -45     0    27    54    27
     0     0     0     0    54     0     0   -54   -54     0
     0     0     0     0    -9     0     0    27     0     0
     0     0    -9     0     9    36     0     0   -27   -27
     0     0     0     0    -9     0     0     0    27     0
     0     0     2     0     0    -9     0     0     0     9
      ] / 2;
 case 4
  V = [
   3   -25   -25    70   140    70   -80  -240  -240   -80    32   128   192   128    32
     0    48     0  -208  -208     0   288   576   288     0  -128  -384  -384  -128     0
     0   -36     0   228    84     0  -384  -432   -48     0   192   384   192     0     0
     0    16     0  -112   -16     0   224    96     0     0  -128  -128     0     0     0
     0    -3     0    22     0     0   -48     0     0     0    32     0     0     0     0
     0     0    48     0  -208  -208     0   288   576   288     0  -128  -384  -384  -128
     0     0     0     0   288     0     0  -672  -672     0     0   384   768   384     0
     0     0     0     0   -96     0     0   480    96     0     0  -384  -384     0     0
     0     0     0     0    16     0     0   -96     0     0     0   128     0     0     0
     0     0   -36     0    84   228     0   -48  -432  -384     0     0   192   384   192
     0     0     0     0   -96     0     0    96   480     0     0     0  -384  -384     0
     0     0     0     0    12     0     0   -48   -48     0     0     0   192     0     0
     0     0    16     0   -16  -112     0     0    96   224     0     0     0  -128  -128
     0     0     0     0    16     0     0     0   -96     0     0     0     0   128     0
     0     0    -3     0     0    22     0     0     0   -48     0     0     0     0    32      
      ]/3;
 otherwise
  % these computations should not be done everytime!
  % for higher pdeg, perhaps also use higher accuracy!!!
  disp(['please hardcode the following matrix into' ...
	' fem_basis_weight_matrix:']);
  % V = inv(P) where P = (p(l1),...p(lm)) powervector of lagrange-nodes
  lagrange_nodes = lagrange_nodes_lcoord(pdeg);
  P = zeros(size(lagrange_nodes,1));
  for i = 1: size(lagrange_nodes,1)
    P(:,i) = power_vector2(lagrange_nodes(i,:),pdeg);
  end;
  V = inv(P); 
end; 
