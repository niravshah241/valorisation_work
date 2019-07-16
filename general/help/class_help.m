function class_help(classname_or_object, verbose)
% function class_help(classname_or_object, verbose)
% prints out class structure and documentation for a given classname or object
%
% The information printed by this function includes
%   - classname,
%   - names of superclasses,
%   - parameters and
%   - methods.
%
%
% Parameters:
%   classname_or_object: Either a string holding the class name or a class
%                        object, for which information shall be printed.
%   verbose            : Integer triggering the amount of information which is
%                        printed.
%                        - '0' - print parameter and method names,
%                        - '1' - print parameter and method names (including
%                                inherited ones),
%                        - '2' - print parameter and method names and their
%                                documentation strings,
%                        - '3' - print parameter and method names and their
%                                documentation strings (including inherited
%                                parameters and methods).

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



  if nargin == 1
    verbose = false;
  end
  if ischar(classname_or_object)
    metainfo = eval(['?',classname_or_object]);
  else
    metainfo = metaclass(classname_or_object);
  end

  cn = metainfo.Name;

  if verbose > 1
    helpstr = @(X) sprintf('\n    %s',strrep(evalc(['help ', X]), sprintf('\n'), sprintf('\n    ')));
  else
    helpstr = @(X) '';
  end

  disp(['CLASS NAME: ', metainfo.Name]);
  disp(' ')
  disp(helpstr(cn));


  inherits = char( cellfun(@(X) [X.Name, ' '], ...
                           metainfo.SuperClasses, ...
                           'UniformOutput', false));
  disp(['Inherits from: ' inherits]);


  if verbose == 0 || verbose == 2
    meths = metainfo.Methods(...
              cellfun(@(X) strcmp(X.DefiningClass.Name,cn)==1, metainfo.Methods, 'UniformOutput',true));
    props = metainfo.Properties(...
              cellfun(@(X) strcmp(X.DefiningClass.Name,cn)==1, metainfo.Properties, 'UniformOutput',true));
  else
    meths = metainfo.Methods;
    props = metainfo.Properties;
  end

  proplist = cellfun(@(X) sprintf('%s%s', ...
                                  X.Name, ...
                                  helpstr([cn, '.', X.Name])), ...
                     props, ...
                     'UniformOutput', false);

  disp(' ')
  disp(' ')
  disp('PROPERTIES');
  disp('==========');
  disp(' ')

  cellfun(@disp, proplist);

  methlist = cellfun(@(X) sprintf('%s%s', ...
                                  X.Name, ...
                                  helpstr([cn, '.', X.Name])), ...
                     meths, ...
                     'UniformOutput', false);

  disp(' ')
  disp(' ')
  disp('METHODS');
  disp('==========');
  disp(' ')

  cellfun(@disp, methlist);
end

