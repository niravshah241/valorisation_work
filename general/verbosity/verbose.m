function r = verbose( level, message, messageId)
% r = verbose( level, message, messageId )
%
% This function displays messages depending on a message-id and/or a level.
% Aditionally you can set/reset the level-filer and add ore remove message-ids
% to the include and exclude lists.
%
% When messages will be displayed?
% --------------------------------
% A Messages only will be diplayed only then if one of the following situation applies.
% 1. The message-id is listet in the include list
% 2. The messages level is lower or equal to the previously set and the message-id don't
%    appear in the exclude list 
% 
% Basic examples:
% ---------------
%  To use the basic operation of verbose, only one or two arguments,
%  the level and the message are necessary.
%
%  three ways to set the verbosity level e.g. to 10
%  1. >  verbose(10);
%  2. > verbose('level',10);
%  3. > oldVerbLev = verbose(10);	% change level to 10
%     > one_of_my_functions();         % ... some code ...
%     > verbose(oldVerbLevel);		% reset level to value before changed
%
%  get the current verbosity level
%  1. > level = verbose();
%  2. > level = verbose('level');
%
%  message output (verbosity level previously set to 10)
%  > verbose(10, 'Hello World!'); % message will be displayed
%  > verbose(11, 'Hello World!'); % ... won't ...
%  > verbose(1,  'Hello World!'); % ... will ...
%
% Extended examples:
% ------------------
% The third argument, the message-id enables you to specify more precisely
% what should be displayed or not.
% 
% add/remove one or more message-ids to the include-/excludelist
% > verbose('addInclude', 'RB:ERROR');                 % add 'RB:ERROR' to the include list
% > verbose('addInclude', {'RB:WARNING', 'RB:INFO'});  % add 'RB:WARNING' and RB:INFO ...
% > verbose('include');                                % return includelist
%                                                      % e.g. {'RB:ERROR', 'RB:WARNING', 'RB:INFO'}
% > verbose('addExclude', {'RB:INFO','RB:UNKNOWN'});   % add 'RB:INFO','RB:UNKNOWN' to excludelist
% > verbose('delInclude', {'RB:ERROR','RB:INFO'});     % remove 'RB:ERROR' and 'RB:INFO'
% > verbose('exclude', {});                            % clear exclude list
%
% print messages 
% > verbose(10);
% > verbose('include', {'RB:SURPRISE1'});
% > verbose('exclude', {'RB:SURPRISE2'});
% > verbose();
% > verbose(1, 'Hello World!', 'RB:INFO');
% > verbose(99, 'Hello What Happened', 'RB:SURPRISE1');   % displayed, cause RB:SURPRISE1 is in includelist
% > verbose(0, 'Surprise for you?');                      % no, you don't get it, it's excluded

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



  persistent VERBOSE;

  % initialize structure VERBOSE
  if isempty( VERBOSE )
    VERBOSE.level = 0;
    VERBOSE.include = {};
    VERBOSE.exclude = {};
    VERBOSE.addInclude = @(x) VERBOSE_addInclude(x);
    VERBOSE.delInclude = @(x) VERBOSE_delInclude(x);
    VERBOSE.addExclude = @(x) VERBOSE_addExclude(x);
    VERBOSE.delExclude = @(x) VERBOSE_delExclude(x);
  end;

  if (nargin==0)
    r = VERBOSE.level;
  elseif (nargin>0) 
    if isnumeric( level )
      if nargin==2
        % filter by level
        if ( level <= VERBOSE.level )
          disp( message );
        end;
        return;
      elseif nargin==3
        [st,i] = dbstack;
        if (length(st)>1)
          [paths,fileName] = fileparts(st(2).file);
          funcName = st(2).name;
          lineNo = num2str(st(2).line);
          if isempty(messageId)
            messageId = sprintf('RB:%s:UNKNOWN',fileName);
          end;
          message = sprintf('%s file=%s func=%s line=%s\n ==> ''%s''\r\n', ...
                             messageId, fileName, funcName, lineNo, message);
        else
          if isempty(messageId)
            messageId = 'RB:UNKNOWN:UNKNOWN';
            message = sprintf('RB:UNKNOWN:UNKNOWN\r\n ==> ''%s''\r\n', message);
          else
            message = sprintf('%s\r\n ==> ''%s''\r\n', messageId, message);
          end;
        end;
        % filter by messageId and level
        if (( ~isempty(VERBOSE.include ) && (  any(strcmp(messageId, VERBOSE.include )))) || ...
            ((level<=VERBOSE.level) && ( isempty(VERBOSE.exclude ) || ( ~any(strcmp(messageId, VERBOSE.exclude )))))) 
          disp( message );
        end;
        return;
      elseif ( nargin == 1 )
        r = VERBOSE.level;
        VERBOSE.level = level;
	return;
      end;
    elseif ischar(level)
      if isfield(VERBOSE, level)
        if isa(VERBOSE.(level), 'function_handle')
          % call function field
          r = VERBOSE.(level)(message);
          return;	
        else
          % set value of field
          r = VERBOSE.(level);
          if nargin==2
            VERBOSE.(level) = message;
          end;
        end;
      else
        verbose( 0, 'first argument is a unknown field.', 'RB:VERBOSE:ARGERR');
      end;
    else
      verbose( 0, 'first argument should be of type numeric or character.', 'RB:VERBOSE:ARGERR');
    end;
  end;

  function r = VERBOSE_addExclude(c)
    r = VERBOSE.exclude;
    VERBOSE.exclude = union(VERBOSE.exclude, c);
  end
    
  function r = VERBOSE_delExclude(c)
    r = VERBOSE.exclude;
    VERBOSE.exclude = setdiff(VERBOSE.exclude, c);
  end

  function r = VERBOSE_addInclude(c)
    r = VERBOSE.include;
    VERBOSE.include = union(VERBOSE.include, c);
  end
    
  function r = VERBOSE_delInclude(c)
    r = VERBOSE.include;
    VERBOSE.include = setdiff(VERBOSE.include, c);
  end

end
%| \docupdate 
