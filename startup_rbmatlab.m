% Script that initializes the RBmatlab-Framework
%
% This script adds all m-files of the RBmatlab package to the Matlab search
% path and checks the environment variables
%  - 'RBMATLABHOME' pointing to the base directory of the RBmatlab
%     installation,
%  - 'RBMATLABTEMP' pointing to a directory where large temporary data can be
%     written to, and the optional
%  - 'RBMATLABRESULT' pointing to a directory where results can be stored.
%  .
% If 'RBMATLABRESULT' is not set, results are written into the directory given
% by 'RBMATLABTEMP'.
%

% get current directory;
fprintf('.__ .__          , .   .  \n')
fprintf('[__)[__)._ _  _.-+-| _.|_ \n')
fprintf('|  \\[__)[ | )(_] | |(_][_)\n')

disp('Starting up RBmatlab in directory:');
p = fileparts(which('startup_rbmatlab'));
disp(p);

% Define the default paths:
paths = {'3rdparty','bin', 'grid', 'general', 'datafunc', 'scripts', 'test', ...
    'demos', 'rbasis', 'discfunc', 'models', 'dune','nirav','prof_haasdonk'};
for i = 1:length(paths)
    addpath(genpath([p, '/', paths{i}]));
end

chdir(p);
clear('p');

% Check if a temporary directory is set and if it exists:
tempdir = getenv('RBMATLABTEMP');
if isempty(tempdir)
  % No RBMATLABTEMP has been provided: Try the default unix /tmp directory:
  tempdir = '/tmp/rbmatlab-temp';
  setenv('RBMATLABTEMP', tempdir);
  status = mkdir(tempdir);
  if ~status      
    error(['Please set an environment-variable RBMATLABTEMP for', ...
         ' temporary data.'])
  end
end
if ~exist(tempdir,'dir')
  % Try to mkdir the dir:
  status = mkdir(tempdir);
  if ~status
    error(['RBMATLABTEMP directory ',tempdir,' does not exist ', ...
        'and cannot be created!']);
  end
end;
disp('Temporary data directory:');
disp(tempdir);

% The result directory:
resultdir = getenv('RBMATLABRESULT');
if isempty(resultdir)
  resultdir = tempdir;
  setenv('RBMATLABRESULT', resultdir);
elseif ~exist(resultdir, 'dir')
  [status,~,~] = mkdir(resultdir);
  if ~status
    error(['RBTMATLABRESULT directory ',resultdir,' does not exist ', ...
        'and cannot be created!']);
  end
end
disp('Result data directory:');
disp(resultdir);

% Home directory:
tdir = [fileparts( which('startup_rbmatlab')),filesep];
homedir = tdir;
setenv('RBMATLABHOME', tdir);
addpath(tdir);
if ~exist(fullfile(tdir,'startup_rbmatlab.m'),'file');
  error(['RBMATLABHOME directory set wrong. No startup_rbmatlab.m' ...
	 ' found.']);
end;
disp('Using the following directory as RBMATLABHOME:');
disp(homedir);

% create cache-directory if not existent:
if ~exist(fullfile(tempdir,'cache'),'dir');
  disp('Creating cache directory')
  success = mkdir(tempdir,'cache');
  if ~success
    error('Error in creating cache subdirectory in temporary dir!');
  end;
end;
clear('tempdir');
clear('tdir');
disp('')
disp('RBmatlab has started successfully!')
