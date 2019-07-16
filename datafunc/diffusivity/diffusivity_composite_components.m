function [C,Cmatrix] = diffusivity_composite_components(glob, params)
% glob is a npoints times 2 matrix
%
% required fields of params:
%       B1, B2: number of composite blocks
%
% generated fields of diffusivity:
%     diffusivity: cell array of vectors with component values
%      { c1,... c_Qc}   with  Qc = B1*B2
%      and c1 = vector with values 0 or 1 of length glob
%
% sanity check: sum_i ci = ones(size(glob,1),1);

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



xB=params.B1;
yB=params.B2;

dx=1/params.B1;
dy=1/params.B2;

nob=xB*yB;

xblocknum = ceil(glob(:,1) /dx);
yblocknum = ceil(glob(:,2) /dy);
blocknum = xblocknum+ params.B1*(yblocknum-1);


C=cell(1,nob);
Cmatrix=[];

for i=1:nob
    %c_zeros=zeros(length(glob),1);
    c_block= blocknum == i;
    %c_zeros(c_block)=1;
    C{1,i}=c_block;
    Cmatrix(:,i)=c_block;
end
