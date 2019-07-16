function diffusivity = diffusivity_composite(glob, params)
%function diffusivity = diffusivity_composite(glob, params)
%
% function computing the diffusivity pointwise evaluation in the point
% sequences indicated by global coordinates in the columns of the matrix glob.
% It returns a piecewise constant diffusion coefficient.
%
% glob is a npoints times 2 matrix
%
% required fields of params:
%       B1, B2: number of composite blocks
%       mu1, mu2,.... mu(B1*B2-1)
%
% generated fields of diffusivity:
%     diffusivity: vector with diffusivity values

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Diffu
EPS = 1e-10;
xB=params.B1;
yB=params.B2;

x_range=0:(1/xB):1;
y_range=0:(1/yB):1;

dx=1/params.B1;
dy=1/params.B2;

nob=xB*yB;

% params.mu_names = {'mu1','mu2',..};
% params.mu1 = 1; 
% params.mu2 = 2;
% params.mu3 = 3;
% params....

% diffu = [1,2,3]
%diffu_value = params.get_mu(params);
diffu_value = [params.get_mu(params), 1];
%diffu_value=params.diffu;
%DIFF_M=flipud(reshape(diffu_value, xB,yB));

params.diffu_value=diffu_value;
%das folgende ohne Schleife:
% % diffu_value = get_mu(model);
% % diffu_value = [mu(:);1];
%dx ist 1/B1 , dy = 1/B2

xblocknum = ceil(glob(:,1) /dx);
yblocknum = ceil(glob(:,2) /dy);
blocknum = xblocknum+ params.B1*(yblocknum-1);
diffusivity = params.diffu_value(blocknum);

diffusivity=diffusivity';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% i=0;
% for i=1:length(com_t)
%     x_check_com_t(i)=com_t{i}(1);
%     y_check_com_t(i)=com_t{i}(2);
%     
%     tmp_x_range_sort=sort([x_range, x_check_com_t(i)]);
%     tmp_y_range_sort=sort([y_range, y_check_com_t(i)]);
%     
%     xBlock(i)=find(tmp_x_range_sort <= (x_check_com_t(i)+EPS) & tmp_x_range_sort >= (x_check_com_t(i)-EPS))-1;
%     yBlock(i)=find(tmp_y_range_sort <= (y_check_com_t(i)+EPS) & tmp_y_range_sort >= (y_check_com_t(i)-EPS))-1;  
%     
%     c_test(i)=DIFF_M(xBlock(i),yBlock(i));
% end























