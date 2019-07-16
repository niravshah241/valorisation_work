% script that initializes the RBmatlab-Framework
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
disp('starting up rbmatlab in directory:');
setpref('Internet','SMTP_Server','localhost');
p = fileparts( which('startup_rbmatlab'));
disp(p);
addpath( fullfile( p ,'bin' ) )
addpath( fullfile( p ,'grid' ) )
%addpath( fullfile( p ,'grid','alu3d' ) )
addpath( fullfile( p ,'grid','common' ) )
addpath( fullfile( p ,'general' ) )
addpath( fullfile( p ,'general','basic' ) )
addpath( fullfile( p ,'general','filecaching' ) )
addpath( fullfile( p ,'general','geometry' ) )
addpath( fullfile( p ,'general','help' ) )
addpath( fullfile( p ,'general','postprocess' ) )
%addpath( fullfile( p ,'general','quadratures' ) )
addpath( fullfile( p ,'general','vecmat' ) )
%addpath( fullfile( p ,'fem' ) )
addpath( fullfile( p ,'datafunc' ) )
addpath( fullfile( p ,'datafunc','auxiliary' ) )
addpath( fullfile( p ,'scripts' ) )
addpath( fullfile( p ,'scripts', 'steps' ) )
addpath( fullfile( p ,'test' ) )
addpath( fullfile( p ,'demos' ) )
addpath( fullfile( p ,'rbasis' ) )
addpath( fullfile( p ,'rbasis','lin_evol' ) )
addpath( fullfile( p ,'rbasis','nonlin_evol' ) )
addpath( fullfile( p ,'rbasis','common' ) )
addpath( fullfile( p ,'rbasis','basisgen' ) )
addpath( fullfile( p ,'rbasis','basisgen','ei' ) )
addpath( fullfile( p ,'rbasis','lin_stat' ) )
addpath( fullfile( p ,'datafiles' ) )
addpath( fullfile( p ,'3rd_party','export_fig' ) )
addpath( fullfile( p ,'discfunc' ) )
addpath( fullfile( p ,'discfunc','common') )
addpath( fullfile( p ,'discfunc','ldg') )
addpath( fullfile( p ,'discfunc','fem') )
addpath( fullfile( p ,'discfunc','fv') )

addpath( fullfile( p, 'rbasis', 'lin_ds') );
addpath( fullfile( p, 'datafunc', 'neumann_values') );
addpath( fullfile( p, 'datafunc', 'init_values') );
addpath( fullfile( p, 'datafunc', 'dirichlet_values') );
addpath( fullfile( p, 'datafunc', 'velocity') );
addpath( fullfile( p, 'datafunc', 'diffusivity') );
addpath( fullfile( p, 'datafunc', 'diffusivity_tensor') );
addpath( fullfile( p, 'datafunc', 'conv_flux') );
addpath( fullfile( p, 'datafunc', 'output_functional') );
addpath( fullfile( p ,'models' ) )
addpath( fullfile( p, 'models', 'advection_output') );
addpath( fullfile( p, 'models', 'common') );
addpath( fullfile( p, 'models', 'convdiff') );
addpath( fullfile( p, 'models', 'dune-rb') );

chdir(p);
clear('p');

tempdir = getenv('RBMATLABTEMP');
if isempty(tempdir) 
  error(['Please set an environment-variable RBMATLABTEMP for', ...
	 ' temporary data'])
end
if ~exist(tempdir,'dir')
  error(['RBMATLABTEMP directory ',tempdir,' does not exist!']);
end;
disp('Using the following directory for large temporary data:');
disp(tempdir);

resultdir = getenv('RBMATLABRESULT');
if isempty(resultdir)
  resultdir = tempdir;
  setenv('RBMATLABRESULT', resultdir);
elseif ~exist(resultdir, 'dir')
  error(['RBMATLABRESULT directory ', resultdir,' does not exist!']);
end
disp('Using the following directory for data files storing results:');
disp(resultdir);

%| \todo What do we need the RBMATLABHOME environment variable for, when the
% addpath commands are executed w.r.t. the current directory?
homedir = getenv('RBMATLABHOME');
if isempty(homedir) 
  error(['Please set an environment-variable RBMATLABHOME pointing', ...
	 ' to the directory RBmatlab'])
end
tdir = [fileparts( which('startup_rbmatlab')),filesep];
if ~exist(fullfile(homedir,'startup_rbmatlab.m'),'file');
  %if ~isequal(tempdir,tdir) => problem with final backlash
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

%disp('clearing filecache for function-calls');
%filecache_clear;
disp('skipped clearing filecache for function-calls!');

%model = convdiff_dune_model('init_model');
%model_data = gen_model_data(model);
%sim_data = detailed_simulation(model, model_data);
%model.mexptr('reset_rb');
%detailed_data.RB = model.mexptr('init_data_basis');
%keyboard;
%detailed_data.RB = model.mexptr('rb_extension_PCA', model.mu, -1, 0.999);
%keyboard;
%reduced_data = gen_reduced_data(model, detailed_data);
%rb_sim_data = rb_simulation(model, reduced_data);

% some change
