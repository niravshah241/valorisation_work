function rv = doxygen(param1, param2)
% function rv = doxygen(param1, param2) is ignored
% Here comes a short description text
%
% After the first empty documentation line, paragraphs of the detailed
% description begin.
%
% @attention If you read the source code of this file, the text might be a bit
% misleading. It is easier and better to read the generated HTML output file.
%
% Here, you have the following formatting possibilities:
% Doxygen commands always begin with an at-character(\@) OR a backslash(\\).
%
% As regular doxygen commands like @verbatim \c @f$, @f$ @f[, @f] @endverbatim
% which produce LaTeX output on HTML websites, look too distracting in matlab
% documentation output, the following shortcut exist: The doxygen filter
% translates
%  - @verbatim 'word' to \c word @endverbatim resulting in the output: 'word',
%  - @verbatim `x` to @f$x@f$ @endverbatim resulting in the output: `x` and 
%  - @verbatim ``x`` to @f[x.@f] @endverbatim resulting in the output: ``x``.
% @note the difference between the two quotation mark for typewriter fonts and
% LaTeX output.
%
% You therefore need to be careful with the use of characters @verbatim ' `
% @endverbatim. If you want to say something about the transposed of a Matrix
% 'A', better do it in a Tex-Environment as `A' * B'` or in a verbatim/code
% environment as
% @code A' * B' @endcode
%
% Words prepended by \\c are written in a \c type-writer font.
% Words prepended by \\b are written in a \b bold font.
% Words prepended by \\em are written in an \em emphasized font.
%
% Blocks starting with @@verbatim or @@code and are ended with @@endverbatim or
% @@endcode are written unformatted in a type-writer font and are not
% interpreted by doxygen.
%
% Example:
% @verbatim
%                /| |\
%               ( |-| )
%                ) " (
%               (>(Y)<)
%                )   (
%               /     \
%              ( (m|m) )  hjw
%            ,-.),___.(,-.\`97
%            \`---\'   \`---\'
% @endverbatim
%
% Paragaphs starting with line ending with a double-colon:
% are started with a bold title line
%
% If you want to reference a function that is elsewhere documented append
% brackets and write @verbatim doxygen() @endverbatim in order to create a link
% for doxygen().
%
% If, however, a double-colon at the end of a line is succeeded by: 
% whitespace characters, like spaces or tabulators the line is not written in a
% bold font.
%
% Listings can be added by prepending lines with a dash(-)
%  - list 1 item which can also include
%   newlines
%  - list 2 item
%    - and they can be nested
%    - subitem 2
%    .
%  - list 3 item
%
% and they are ended by an empty documentation line.
%
% Enumerations can be added by prepending lines with a dash and hash (-#)
%  -# first item
%  -# second item
%  -# third item
%
% Lines beginning with the words "Parameters" or "Return values" start a block
% of argument respectively return argument descriptions.
%
% Parameters:
%  param1: first parameter
%
% Return values:
%  rv: return value
%
% A line beginning with the words "Required fields of", "optional fields of" or
% "generated fields of" start a block for descriptions for fields used by the
% parameters or generated for the return values.
%
% Required fields of param1:
%  test: Description for required field param1.test
%
% Optional fields of param2:
%  test2: Description for optional field param2.test2
%
% Generated fields of rv:
%  RB: Description for generated field rv.RB
%
%

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


% After the first non-comment line the function body begins:

%| Comment blocks starting with %| are interpreted as Doxygen documentation
% blocks and can include doxygen commands like

%| \todo There needs to be done something in this file

% fields of parameters that are used in the function body are added to the
% required fileds list automatically, if they are not documented yet.
param1.auto_added;

param2.auto_added;

% fields of return values that are assigned somewhere in the function body are
% also added automatically to the list of generated fields
rv.auto_added  = 1;
rv.sub.auto_added = 2;

param1.sub.auto_added;

