function plot_basisgen_results(fns, model)
%function plot_basisgen_results(fns[, model])
%
% function plotting of basis generation results generated
% by basisgen_main. Fns is a cell array with names of the
% basisgen-files assumed to exist in $(RBMATLABTEMP)/basisgen
% Some fields in model can be used to select the produced plots.
% by default all plots are generated
%
% required fields of model:
%         plot_train_times: flag indicating, whether the train times are
%                     to be compared  
%         plot_train_estimators: flag indicating, whether the train
%                     error estimators are to be compared  
%         plot_test_estimators: flag indicating, whether the test
%                     error estimators are to be compared  
%         plot_test_errors: flag indicating, whether the test
%                     error estimators are to be compared  
%         plot_legends : cell array of legends to be used
%         plot_linestyles : cell array of linestyles to be used
%         plot_linecolors : cell array of linecolors
%         plot_linewidths : vector of linewidths

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


% Bernard Haasdonk 7.6.2007

if nargin<2
  model = [];
end;

if ~isfield(model,'plot_train_times')
  model.plot_train_times = 1;
end;

if ~isfield(model,'plot_train_estimators')
  model.plot_train_estimators = 1;
end;

if ~isfield(model,'plot_test_estimators')
  model.plot_test_estimators = 1;
end;

if ~isfield(model,'plot_test_errors')
  model.plot_test_errors = 1;
end;

if ~isfield(model,'fontsize')
  model.fontsize = 10;
end;

if ~isfield(model,'plot_legends_location')
  legends_location = 'NorthEast';
else
  legends_location = model.plot_legends_location;
end;

%fns = {'RB_random_fixed_50_tested',...
%       'RB_random_fixed_100_tested',...
%       'RB_uniform_fixed_49_tested',...
%       'RB_uniform_fixed_100_tested',...
%       'RB_uniform_refined_3_3_tested',...
%       'RB_adaptive_refined_3_3_tested'};
%fns = {       'RB_random_fixed_100_tested',...
%              'RB_uniform_fixed_100_tested',...
%       'RB_uniform_refined_3_3_tested',...
%       'RB_adaptive_refined_3_3_tested'};
%fns = {'RB_uniform_fixed_100_tested',...
%       'RB_uniform_refined_3_3_tested',...
%       'RB_adaptive_refined_3_3_tested'};
%fns = {'RB_uniform_fixed_100',...
%       'RB_adaptive_refined_2_2'};

%fns= {'RB_uniform_fixed_16_tested', ...
%		    'RB_uniform_fixed_64_tested', ...
%		    'RB_uniform_fixed_256_tested', ...
%		    'RB_uniform_fixed_1024_tested'};

if isempty(fns)  
  fns = { ...
      'RB_random_fixed_16_tested', ...
      'RB_random_fixed_64_tested', ...
      'RB_random_fixed_256_tested', ...
      'RB_random_fixed_1024_tested'};
end;
%fns = {'RB_adaptive_refined_2_2_r1_tested', ...
%		    'RB_adaptive_refined_2_2_r2_tested', ...
%		    'RB_adaptive_refined_2_2_r5_tested', ...
%		    'RB_adaptive_refined_2_2_r10_tested'};

%fns = {'RB_uniform_fixed_49',...
%       'RB_uniform_fixed_longt_49'};
%compute_RB_methods = {'RB_random_fixed_100'};

%fns = {'RB_adaptive_refined_tested',...
%       'RB_uniform_fixed_tested','RB_uniform_refined_tested'};

data = {};
for fnum = 1:length(fns);
  data{fnum} = load(fullfile(getenv('PDEMATLABTEMP'),'basisgen',fns{fnum})); 
end;



%         plot_legends : cell array of legends to be used
%         plot_linestyles : cell array of linestyles to be used
%         plot_linecolors : cell array of linecolors
%         plot_linewidths : vector of linewidths

if ~isfield(model,'plot_legends');
  % generate legends (without final '_tested'):
  lfns = fns;
  for fnum = 1:length(fns);
    fn = fns{fnum};
    i = findstr(fn,'_');
    for j = length(i):-1:1
      fn = [fn(1:(i(j)-1)),'\_',fn(i(j)+1:end)];
    end;
    i = findstr(fn,'\_tested');
    if ~isempty(i)
      fn = fn(1:(i(1)-1));
    end;
    lfns{fnum} = fn;
  end;
else
  lfns = model.plot_legends;
end;


if model.plot_train_times
  
  % plot of computation times
  times = [];
  for fnum = 1:length(fns);
    l = length(data{fnum}.detailed_data.RB_info.toc_value_sequence);
    if size(times,2)> l
      times = [times;data{fnum}.detailed_data.RB_info.toc_value_sequence, ...
	       NaN * ones(1,size(times,2)-l)];
    else
      times = [times, NaN * ones(size(times,1),l-size(times,2));...
	       data{fnum}.detailed_data.RB_info.toc_value_sequence];    
    end;
  end;
  figure; 
  p = plot(times');
  for i = 1:length(p)
    if isfield(model,'plot_linecolors')
      set(p(i),'color',model.plot_linecolors{i});
    end;
    if isfield(model,'plot_linewidths')
      set(p(i),'Linewidth',model.plot_linewidths(i));
    end;
    if isfield(model,'plot_linestyles')
      set(p(i),'Linestyle',model.plot_linestyles{i});
    end;
  end;
  
  title('computation time')
  xlabel('num basis functions N');
  ylabel('computation time [s]')
  legend(lfns,'location',legends_location);
end;

if model.plot_train_estimators
  % plot of train-indicator decrease
  estimators = [];
  for fnum = 1:length(fns);
    l = length(data{fnum}.detailed_data.RB_info.max_err_sequence);
    if size(estimators,2)> l
      estimators = [estimators;data{fnum}.detailed_data.RB_info.max_err_sequence, ...
		    NaN * ones(1,size(estimators,2)-l)];
    else
      estimators = [estimators, NaN * ones(size(estimators,1),...
					   l-size(estimators,2));...
		    data{fnum}.detailed_data.RB_info.max_err_sequence];    
    end;
  end;
  figure, p = plot(estimators');
  for i = 1:length(p)
    if isfield(model,'plot_linecolors')
      set(p(i),'color',model.plot_linecolors{i});
    end;
    if isfield(model,'plot_linewidths')
      set(p(i),'Linewidth',model.plot_linewidths(i));
    end;
    if isfield(model,'plot_linestyles')
      set(p(i),'Linestyle',model.plot_linestyles{i});
    end;
  end;

  
  title('train error-indicator value decrease')
  xlabel('num basis functions N');
  ylabel('maximum indicator values')
  set(gca,'Yscale','log');
  legend(lfns,'location',legends_location);
end;

max_test_error_available = zeros(length(fns));
max_test_estimator_available = zeros(length(fns));
min_test_estimator_available = zeros(length(fns));
min_test_error_available = zeros(length(fns));
for fnum = 1:length(fns);
  if isfield(data{fnum}.detailed_data.RB_info, ...
	     'max_test_error_sequence')
    max_test_error_available(fnum) = 1;
  end;
  if isfield(data{fnum}.detailed_data.RB_info, ...
	     'max_test_estimator_sequence')
    max_test_estimator_available(fnum) = 1;
  end;
  if isfield(data{fnum}.detailed_data.RB_info, ...
	     'min_test_estimator_sequence')
    min_test_estimator_available(fnum) = 1;
  end;
  if isfield(data{fnum}.detailed_data.RB_info, ...
	     'min_test_error_sequence')
    min_test_error_available(fnum) = 1;
  end;
end;

%--------------------------------------------------------------------------

% plot of max test-estimator decrease
if model.plot_test_estimators
  estimators = [];
  estimatorsX = [];
  estimatorsY = [];
  fi = find(max_test_estimator_available);
  for fnum = 1:length(fi);
    max_seq = data{fi(fnum)}.detailed_data.RB_info.max_test_estimator_sequence;
%    l = length(max_seq);
%    if size(estimators,2)> l
%      estimators = ...
%	  [estimators;...
%	   max_seq(:)',...
%	   NaN * ones(1,size(estimators,2)-l)];
%    else
%      estimators = [estimators, ...
%		    NaN * ones(size(estimators,1),...
%			       l-size(estimators,2));...
%		    max_seq(:)'];
%    end;
    
    max_seqX = find(~isnan(max_seq));
    max_seqY = max_seq(max_seqX);
    l = length(max_seqX);
    if size(estimatorsX,2)> l
      estimatorsX = ...
	  [estimatorsX;...
	   max_seqX(:)',...
	   NaN * ones(1,size(estimatorsX,2)-l)];
      estimatorsY = ...
	  [estimatorsY;...
	   max_seqY(:)',...
	   NaN * ones(1,size(estimatorsY,2)-l)];
    else
      estimatorsX = [estimatorsX, ...
		     NaN * ones(size(estimatorsX,1),...
				l-size(estimatorsX,2));...
		     max_seqX(:)'];
      estimatorsY = [estimatorsY, ...
		     NaN * ones(size(estimatorsY,1),...
				l-size(estimatorsY,2));...
		     max_seqY(:)'];
    end;    
    
  end;
  if ~isempty(fi)
    figure;
    %    if isempty(isnan(estimators))
    %      p = plot(estimators');
    %    else
    %      p = plot(estimators','x');
    %    end;
    p = plot(estimatorsX',estimatorsY','x');
    for i = 1:length(p)
      if isfield(model,'plot_linecolors')
	set(p(i),'color',model.plot_linecolors{fi(i)});
      end;
      if isfield(model,'plot_linewidths')
	set(p(i),'Linewidth',model.plot_linewidths(fi(i)));
      end;
      if isfield(model,'plot_linestyles')
	set(p(i),'Linestyle',model.plot_linestyles{fi(i)});
      end;
    end;
    
    title('test estimator values decrease','Fontsize',model.fontsize)
    xlabel('num basis functions N','Fontsize',model.fontsize);
    ylabel('maximum test estimator values Delta_N','Fontsize',model.fontsize)
    set(gca,'Yscale','log');
    legend(lfns(fi),'location',legends_location,'Fontsize',model.fontsize);
  end;
end;


% plot of min test-estimator decrease
if model.plot_test_estimators
  estimators = [];
  fi = find(min_test_estimator_available);
  for fnum = 1:length(fi);
    l = length(data{fi(fnum)}.detailed_data.RB_info.min_test_estimator_sequence);
    if size(estimators,2)> l
      estimators = ...
	  [estimators;...
	   data{fi(fnum)}.detailed_data.RB_info.min_test_estimator_sequence(:)', ...
	   NaN * ones(1,size(estimators,2)-l)];
    else
      estimators = [estimators, ...
		    NaN * ones(size(estimators,1),...
			       l-size(estimators,2));...
		    data{fi(fnum)}.detailed_data.RB_info.min_test_estimator_sequence(:)'];    
    end;
  end;
  if ~isempty(fi)
    figure;
    if isempty(isnan(estimators))
      p = plot(estimators');
    else
      p = plot(estimators','x');
    end;
    for i = 1:length(p)
      if isfield(model,'plot_linecolors')
	set(p(i),'color',model.plot_linecolors{fi(i)});
      end;
      if isfield(model,'plot_linewidths')
	set(p(i),'Linewidth',model.plot_linewidths(fi(i)));
      end;
      if isfield(model,'plot_linestyles')
	set(p(i),'Linestyle',model.plot_linestyles{fi(i)});
      end;
    end;
    title('test estimator values decrease')
    xlabel('num basis functions N');
    ylabel('minimum test estimator values Delta_N')
    set(gca,'Yscale','log');
    legend(lfns(fi),'location',legends_location);
  end;
end;


% plot of max/min test-estimator ratio
if model.plot_test_estimators
  estimators = [];
  
  estimatorsX = [];
  estimatorsY = [];

  for fnum = 1:length(fi);
    max_seq = data{fi(fnum)}.detailed_data.RB_info.max_test_estimator_sequence;
    min_seq = data{fi(fnum)}.detailed_data.RB_info.min_test_estimator_sequence;
    seqX = find(~isnan(max_seq));
    seqY = max_seq(seqX)./min_seq(seqX);
    l = length(seqX);
    if size(estimatorsX,2)> l
      estimatorsX = ...
	  [estimatorsX;...
	   seqX(:)',...
	   NaN * ones(1,size(estimatorsX,2)-l)];
      estimatorsY = ...
	  [estimatorsY;...
	   seqY(:)',...
	   NaN * ones(1,size(estimatorsY,2)-l)];
    else
      estimatorsX = [estimatorsX, ...
		     NaN * ones(size(estimatorsX,1),...
				l-size(estimatorsX,2));...
		     seqX(:)'];
      estimatorsY = [estimatorsY, ...
		     NaN * ones(size(estimatorsY,1),...
				l-size(estimatorsY,2));...
		     seqY(:)'];
    end;    

  end;
  if ~isempty(fi)
    figure;
    %    if isempty(isnan(estimators))
    %      p = plot(estimators');
    %    else
    %      p = plot(estimators','x');
    %    end;
    
    p = plot(estimatorsX',estimatorsY','x');
    
    for i = 1:length(p)
      if isfield(model,'plot_linecolors')
	set(p(i),'color',model.plot_linecolors{fi(i)});
      end;
      if isfield(model,'plot_linewidths')
	set(p(i),'Linewidth',model.plot_linewidths(fi(i)));
      end;
      if isfield(model,'plot_linestyles')
	set(p(i),'Linestyle',model.plot_linestyles{fi(i)});
      end;
    end;
    title('max/min test estimator ratio','Fontsize',model.fontsize)
    xlabel('num basis functions N','Fontsize',model.fontsize);
    ylabel('max/min test estimator ratio','Fontsize',model.fontsize)
    set(gca,'Yscale','log');
    legend(lfns(fi),'location',legends_location,'Fontsize',model.fontsize);
  end;
end;


if isfield(model,'plot_max_test_vs_time')
  
  % plot of max test-estimator versus training time 
  estimators = [];  
  estimatorsX = [];
  estimatorsY = [];
  
  for fnum = 1:length(fi);
    max_seq = data{fi(fnum)}.detailed_data.RB_info.max_test_estimator_sequence;
    seqN = find(~isnan(max_seq));
    seqX = data{fi(fnum)}.detailed_data.RB_info.toc_value_sequence(seqN);
    seqY = max_seq(seqN);
    l = length(seqX);
    if size(estimatorsX,2)> l
      estimatorsX = ...
	  [estimatorsX;...
	   seqX(:)',...
	   NaN * ones(1,size(estimatorsX,2)-l)];
      estimatorsY = ...
	  [estimatorsY;...
	   seqY(:)',...
	   NaN * ones(1,size(estimatorsY,2)-l)];
    else
      estimatorsX = [estimatorsX, ...
		     NaN * ones(size(estimatorsX,1),...
				l-size(estimatorsX,2));...
		     seqX(:)'];
      estimatorsY = [estimatorsY, ...
		     NaN * ones(size(estimatorsY,1),...
				l-size(estimatorsY,2));...
		     seqY(:)'];
    end;    
    
  end;
  if ~isempty(fi)
    figure;
    %    if isempty(isnan(estimators))
    %      p = plot(estimators');
    %    else
    %      p = plot(estimators','x');
    %    end;
    
    p = plot(estimatorsX',estimatorsY','x');
    
    for i = 1:length(p)
      if isfield(model,'plot_linecolors')
	set(p(i),'color',model.plot_linecolors{fi(i)});
      end;
      if isfield(model,'plot_linewidths')
	set(p(i),'Linewidth',model.plot_linewidths(fi(i)));
      end;
      if isfield(model,'plot_linestyles')
	set(p(i),'Linestyle',model.plot_linestyles{fi(i)});
      end;
    end;
    title('max test estimator over training time','Fontsize',model.fontsize)
    xlabel('training time','Fontsize',model.fontsize);
    ylabel('max test estimator value','Fontsize',model.fontsize)
    set(gca,'Yscale','log');
    legend(lfns(fi),'location',legends_location,'Fontsize',model.fontsize);
  end;
end;

% plot of max test-error decrease
if model.plot_test_errors
  estimators = [];
  fi = find(max_test_error_available);
  for fnum = 1:length(fi);
    l = length(data{fi(fnum)}.detailed_data.RB_info.max_test_error_sequence);
    if size(estimators,2)> l
      estimators = [estimators;...
		    data{fi(fnum)}.detailed_data.RB_info.max_test_error_sequence(:)',...
		    NaN * ones(1,size(estimators,2)-l)];
    else
      estimators = ...
	  [estimators, NaN * ones(size(estimators,1),l-size(estimators,2));...
	   data{fi(fnum)}.detailed_data.RB_info.max_test_error_sequence(:)'];    
    end;
  end;
  if ~isempty(fi)
    figure;
    if isempty(isnan(estimators))
      p = plot(estimators');
    else
      p = plot(estimators','x');
    end;
    for i = 1:length(p)
      if isfield(model,'plot_linecolors')
	set(p(i),'color',model.plot_linecolors{fi(i)});
      end;
      if isfield(model,'plot_linewidths')
	set(p(i),'Linewidth',model.plot_linewidths(fi(i)));
      end;
      if isfield(model,'plot_linestyles')
	set(p(i),'Linestyle',model.plot_linestyles{fi(i)});
      end;
    end;
    title('test error values decrease')
    xlabel('num basis functions N');
    ylabel('maximum test error values')
    set(gca,'Yscale','log');
    legend(lfns(fi),'location',legends_location);
  end;
end;



% plot of min test-error decrease
if model.plot_test_errors
  estimators = [];
  fi = find(min_test_error_available);
  for fnum = 1:length(fi);
    l = length(data{fi(fnum)}.detailed_data.RB_info.min_test_error_sequence);
    if size(estimators,2)> l
      estimators = [estimators;...
		    data{fi(fnum)}.detailed_data.RB_info.min_test_error_sequence(:)',...
		    NaN * ones(1,size(estimators,2)-l)];
    else
      estimators = ...
	  [estimators, NaN * ones(size(estimators,1),l-size(estimators,2));...
	   data{fi(fnum)}.detailed_data.RB_info.min_test_error_sequence(:)'];    
    end;
  end;
  if ~isempty(fi)
    figure;
    if isempty(isnan(estimators))
      p = plot(estimators');
    else
      p = plot(estimators','x');
    end;
    for i = 1:length(p)
      if isfield(model,'plot_linecolors')
	set(p(i),'color',model.plot_linecolors{fi(i)});
      end;
      if isfield(model,'plot_linewidths')
	set(p(i),'Linewidth',model.plot_linewidths(fi(i)));
      end;
      if isfield(model,'plot_linestyles')
	set(p(i),'Linestyle',model.plot_linestyles{fi(i)});
      end;
    end;
    title('test error values decrease')
    xlabel('num basis functions N');
    ylabel('minimum test error values')
    set(gca,'Yscale','log');
    legend(lfns(fi),'location',legends_location);
  end;
end;





