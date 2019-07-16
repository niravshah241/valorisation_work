function p = plot_mu_frequency(mu_values,params)
%function p = plot_mu_frequency(mu_values,params)
%
% function plotting the mu_values (columnwise 2D/3D vectors) as a
% point distribution, the colors and size indicating the frequency
% of the points
% 
% required fields of params:
%    mu_names: cell array of the names of the components of mu_values

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


% determine frequency vectors in a list

% marker sizes minimum and maximum
mmin = 2;
mmax = 20;

% eliminate possible NaN entries in list:

while ~isempty(find(isnan(mu_values(:,1))))
  mu_values = [mu_values(:,2:end)]; 
end;

if (length(params.mu_names) == 2)
  mu_values = [mu_values; zeros(1,size(mu_values,2))];
  params.mu_names = [params.mu_names, {''}];
end;

[mu_single, mu_frequency] = count_column_frequency(mu_values);

%l = plot3(mu_values(1,:),mu_values(2,:),mu_values(3,:),'.');

cmax = max(mu_frequency);

%cm = jet(cmax);
cm = gray(cmax);

plot_args = {};
pfreq = [];
figure;
for i=1:cmax
  j = find(mu_frequency==i);
  if ~isempty(j)
    plot_args = [plot_args, ...
		 {mu_single(1,j)},{mu_single(2,j)},{mu_single(3, ...
						  j)},{'.'}];
    pfreq = [pfreq, i];
  end;
end;
% first plot all black
%p = plot3(plot_args{:}), hold on;
%for i=1:length(p)
%  wi = (mmax*((pfreq(i)-1)/(cmax-1)) + mmin*(1-(pfreq(i)-1)/(cmax-1)))*1.1;
%  set(p(i),'Markersize',wi,'Color',cm(1,:));
%end;
% then plot all with gray shades
p = plot3(plot_args{:});
for i=1:length(p)
  wi = mmax*((pfreq(i)-1)/(cmax-1)) + mmin*(1-(pfreq(i)-1)/(cmax-1));
  set(p(i),'Markersize',wi,'MarkerFaceColor',cm(pfreq(i),:),'Marker','o',...
	   'MarkerEdgeColor',[0 0 0]);
end;
set(gca,'Box','on');

c = colorbar;
colormap(cm);
set(c,'Clim',[0.5,cmax+0.5]);
c = colorbar;

xlabel(params.mu_names{1});
ylabel(params.mu_names{2});
zlabel(params.mu_names{3});

% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
