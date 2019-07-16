function [x,y]=decomposed_grid(bs,s)
%DECOMPOSED_GRID	Gives geometry data for the decomposed_grid PDE model.
%
%   NE=DECOMPOSED_GRID gives the number of boundary segments
%
%   D=DECOMPOSED_GRID(BS) gives a matrix with one column for each boundary segment
%   specified in BS.
%   Row 1 contains the start parameter value.
%   Row 2 contains the end parameter value.
%   Row 3 contains the number of the left-hand regions.
%   Row 4 contains the number of the right-hand regions.
%
%   [X,Y]=DECOMPOSED_GRID(BS,S) gives coordinates of boundary points. BS specifies the
%   boundary segments and S the corresponding parameter values. BS may be
%   a scalar.

nbs=16;

if nargin==0,
  x=nbs; % number of boundary segments
  return
end

d=[
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 % start parameter value
  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 % end parameter value
  1 4 0 0 4 3 0 0 2 3 1 2 0 0 0 0 % left hand region
  0 0 1 2 0 0 2 3 1 4 4 3 1 4 3 2 % right hand region
];

bs1=bs(:)';

if find(bs1<1 | bs1>nbs),
  error(message('pde:wgeom:NonExistBoundSeg'))
end

if nargin==1,
  x=d(:,bs1);
  return
end

x=zeros(size(s));
y=zeros(size(s));
[m,n]=size(bs);
if m==1 & n==1,
  bs=bs*ones(size(s)); % expand bs
elseif m~=size(s,1) | n~=size(s,2),
  error(message('pde:wgeom:BsSizeError'));
end

if ~isempty(s),

% boundary segment 1
ii=find(bs==1);
if length(ii)
x(ii)=(0.5-(0))*(s(ii)-d(1,1))/(d(2,1)-d(1,1))+(0);
y(ii)=(0-(0))*(s(ii)-d(1,1))/(d(2,1)-d(1,1))+(0);
end

% boundary segment 2
ii=find(bs==2);
if length(ii)
x(ii)=(1-(0.5))*(s(ii)-d(1,2))/(d(2,2)-d(1,2))+(0.5);
y(ii)=(0-(0))*(s(ii)-d(1,2))/(d(2,2)-d(1,2))+(0);
end

% boundary segment 3
ii=find(bs==3);
if length(ii)
x(ii)=(0-(0))*(s(ii)-d(1,3))/(d(2,3)-d(1,3))+(0);
y(ii)=(0.5-(0))*(s(ii)-d(1,3))/(d(2,3)-d(1,3))+(0);
end

% boundary segment 4
ii=find(bs==4);
if length(ii)
x(ii)=(0-(0))*(s(ii)-d(1,4))/(d(2,4)-d(1,4))+(0);
y(ii)=(1-(0.5))*(s(ii)-d(1,4))/(d(2,4)-d(1,4))+(0.5);
end

% boundary segment 5
ii=find(bs==5);
if length(ii)
x(ii)=(1-(1))*(s(ii)-d(1,5))/(d(2,5)-d(1,5))+(1);
y(ii)=(0.5-(0))*(s(ii)-d(1,5))/(d(2,5)-d(1,5))+(0);
end

% boundary segment 6
ii=find(bs==6);
if length(ii)
x(ii)=(1-(1))*(s(ii)-d(1,6))/(d(2,6)-d(1,6))+(1);
y(ii)=(1-(0.5))*(s(ii)-d(1,6))/(d(2,6)-d(1,6))+(0.5);
end

% boundary segment 7
ii=find(bs==7);
if length(ii)
x(ii)=(0.5-(0))*(s(ii)-d(1,7))/(d(2,7)-d(1,7))+(0);
y(ii)=(1-(1))*(s(ii)-d(1,7))/(d(2,7)-d(1,7))+(1);
end

% boundary segment 8
ii=find(bs==8);
if length(ii)
x(ii)=(1-(0.5))*(s(ii)-d(1,8))/(d(2,8)-d(1,8))+(0.5);
y(ii)=(1-(1))*(s(ii)-d(1,8))/(d(2,8)-d(1,8))+(1);
end

% boundary segment 9
ii=find(bs==9);
if length(ii)
x(ii)=(0.29999999999999999-(0))*(s(ii)-d(1,9))/(d(2,9)-d(1,9))+(0);
y(ii)=(0.5-(0.5))*(s(ii)-d(1,9))/(d(2,9)-d(1,9))+(0.5);
end

% boundary segment 10
ii=find(bs==10);
if length(ii)
x(ii)=(1-(0.69999999999999996))*(s(ii)-d(1,10))/(d(2,10)-d(1,10))+(0.69999999999999996);
y(ii)=(0.5-(0.5))*(s(ii)-d(1,10))/(d(2,10)-d(1,10))+(0.5);
end

% boundary segment 11
ii=find(bs==11);
if length(ii)
x(ii)=(0.5-(0.5))*(s(ii)-d(1,11))/(d(2,11)-d(1,11))+(0.5);
y(ii)=(0.29999999999999999-(0))*(s(ii)-d(1,11))/(d(2,11)-d(1,11))+(0);
end

% boundary segment 12
ii=find(bs==12);
if length(ii)
x(ii)=(0.5-(0.5))*(s(ii)-d(1,12))/(d(2,12)-d(1,12))+(0.5);
y(ii)=(1-(0.69999999999999996))*(s(ii)-d(1,12))/(d(2,12)-d(1,12))+(0.69999999999999996);
end

% boundary segment 13
ii=find(bs==13);
if length(ii)
x(ii)=0.20000000000000001*cos(1.5707963267948966*s(ii)+(-3.1415926535897931))+(0.5);
y(ii)=0.20000000000000001*sin(1.5707963267948966*s(ii)+(-3.1415926535897931))+(0.5);
end

% boundary segment 14
ii=find(bs==14);
if length(ii)
x(ii)=0.20000000000000001*cos(1.5707963267948966*s(ii)+(-1.5707963267948966))+(0.5);
y(ii)=0.20000000000000001*sin(1.5707963267948966*s(ii)+(-1.5707963267948966))+(0.5);
end

% boundary segment 15
ii=find(bs==15);
if length(ii)
x(ii)=0.20000000000000001*cos(1.5707963267948966*s(ii)+(0))+(0.5);
y(ii)=0.20000000000000001*sin(1.5707963267948966*s(ii)+(0))+(0.5);
end

% boundary segment 16
ii=find(bs==16);
if length(ii)
x(ii)=0.20000000000000001*cos(1.5707963267948966*s(ii)+(1.5707963267948966))+(0.5);
y(ii)=0.20000000000000001*sin(1.5707963267948966*s(ii)+(1.5707963267948966))+(0.5);
end

end
