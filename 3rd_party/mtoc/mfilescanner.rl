#include "mfilescanner.h"

#include <cassert>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <iterator>
extern "C" {
#include <unistd.h>
}

using std::cerr;
using std::cout;
using std::cin;
using std::endl;
using std::string;
using std::vector;
using std::list;
using std::copy;
using std::map;
using std::set;
using std::istream;
using std::ostream_iterator;
using std::ostringstream;


%%{
  machine MFileScanner;
  write data;

  # end of file character
  EOF = 0;

  # any character other than end of file
  default = ^0;

  # end of line character
  EOL = ('\r'? . '\n') @{line++;};

  # scanner for comment blocks
  in_comment_block :=
  (
   # comment line begins with a percent sign
   '%'
     @{ tmp_p = p+1; cout << " *"; }
   # and then some default characters
   . (default - '\n')* . EOL
     @{ cout.write(tmp_p, p - tmp_p+1); }
  )*
  $!{
    cout << "*/\n";
    if(is_getter_ || is_setter_)
    {
      cout << "/* ";
    }
    fhold;
    fret;
  };

  action end_doxy_block
  {
    if(!docline)
    {
      p = ts-1;
      if(is_class_)
      {
        if(class_part_ == Header)
        {
          end_of_class_doc();
          fgoto classbody;
        } else if(class_part_ == Method || class_part_ == AtMethod)
        {
          print_function_synopsis();
          fgoto funcbody;
        }
        else if(class_part_ == MethodDeclaration)
        {
          fgoto funcdef;
        }
        else if(class_part_ == Property)
        {
          fgoto propertybody;
        }
        else if(class_part_ == InClassComment)
        {
          class_part_ = Method;
          fgoto methods;
        }
        else
        {
          cerr << "missing class part handling for class part: " << ClassPartNames[class_part_] << endl;
        }
      }
      else
      {
        print_function_synopsis();
        fgoto funcbody;
      }
    }
  }

  # executed when end of file is reached
  action end_of_file
  {
    end_function();
    for(  list<string>::iterator it = namespaces_.begin();
          it != namespaces_.end(); ++it)
    {
      cout << "};\n";
    }
  }

  # executed when we reached a comment block
  action in_c_block
  {
    assert(p >= tmp_p-1);
    cout.write(tmp_p, p-tmp_p+1);
    fcall in_comment_block;
  }

  action echo { cout << fc; }

  action st_tok { tmp_p = p; }

  action echo_tok {
    assert (p >= tmp_p);
    cout.write(tmp_p, p - tmp_p);
  }

  action string_tok {
    assert ( p >= tmp_p );
    tmp_string.assign(tmp_p, p-tmp_p);
  }

  # common definitions {{{2

  # comment in function body that might also be added to the doxygen block for
  # the function description
  is_doxy_comment =
    (
    # if percent character is followed by a bar we make the comment a doxygen
    # comment
     '|' @{ if(is_getter_ || is_setter_)
            {
              cout << "*/";
            }
            cout << "/**"; tmp_p = p+1;
          }
     . (default - '\n')* . EOL
     |
    # else: a regular comment
     ( (default - '|')
       @{
         if(is_getter_ || is_setter_)
         {
           cout << "*/";
         }
         cout << "/* ";
         tmp_p = p;
         } )
     . (default - '\n')* . EOL
    );

  # comment block in function body
  comment_block = [ \t]* . '%' . is_doxy_comment;

  # an empty line
  empty_line = [\t ]* . EOL;

  # documentation line begin
  doc_begin = [\t ]* . '%' @{ tmp_p = p + 1; };

  # swallow a comment line till the end of the line (make it a c comment)
  garble_comment_line =
    ( (default - [\r\n])* . EOL )
      @{
        if(is_getter_ || is_setter_)
        {
          cout << "*/";
        }
        cout << "/* ";
        assert( p >= tmp_p );
        cout.write(tmp_p, p - tmp_p) << "*/\n";
        if(is_getter_ || is_setter_)
        {
          cout << "/* ";
        }
      };
  garble_comment_line_wo_eol =
    (default - [\r\n])*;

  # white space or comment
  WSOC =
    ( [ \t]+
      | ('%' @{ tmp_p = p+1; } . garble_comment_line)
      | ('...'.[ \t]*.EOL)
    );

  # matlab identifier
  IDENTEND = [A-Za-z0-9_];
  IDENT = [A-Za-z_]IDENTEND**;


  # matlab identifier with .
  IDENT_W_DOT = [A-Za-z_][A-Za-z0-9_.]**;

  # default arguments in function declarations
  default_arg = [^,)]** @echo;

  #}}}2

  # parameter list for functions {{{2
  paramlist =
    (
     (WSOC | [,\n]
#       @{if(*p=='\n' || paramlist_.size() != 1 || paramlist_[0] != string("this" )) {
#         //buffer_.append(std::string(*p));
#         } }
     | ( '=' . default_arg ) )+
     |
     # matlab identifier (parameter)
     (IDENT)
       >st_tok
       %{
         assert(p >= tmp_p);
         string s(tmp_p, p - tmp_p);
         bool addBlock = true;
         // do not print this pointer
         if( is_class_ && ( ( class_part_ == Method
                              && cfuncname_ != classname_
                              && !methodparams_.statical
                            )
                            || class_part_ == AtMethod
                            || class_part_ == MethodDeclaration
                          )
           )
         {
            if(paramlist_.empty())
            {
              addBlock = false;
              paramlist_.push_back(string("this"));
            }
            else if(paramlist_.size() == 1 && paramlist_[0] == string("this"))
              paramlist_.clear();
         }

         if(addBlock) {

#ifdef DEBUG
{
  ostringstream oss;
  oss << "found parameter: " << s;
  debug_output(oss.str(), p);
}
#endif
           // add an empty docu block for parameter \a s
           param_list_[s] = DocuBlock();
#ifdef DEBUG
{
  ostringstream oss;
  oss << "in paramlist: add to paramlist: " << s;
  debug_output(oss.str(), p);
}
#endif
           paramlist_.push_back(s);
         }
       }
    )**;

  # return parameter list for functions
  lparamlist =
    ( (WSOC | [\n])+
      | ','
      # matlab identifier (return value)
      | ( IDENT > st_tok
          %{
            assert(p >= tmp_p);
            string s(tmp_p, p - tmp_p);
            returnlist_.push_back(s);
            // add an empty docu block for return value \a s
            return_list_[s] = DocuBlock();
          }
        )
    )**;

  # return parameter or return parameter list
  lparams =
    (
      (
        (
         # matlab identifier
         IDENT
           >st_tok
           %{
             assert(p >= tmp_p);
             string s(tmp_p, p - tmp_p);
             returnlist_.push_back(s);
             // add an empty docu block for single return value \a s
             return_list_[s] = DocuBlock();
#ifdef DEBUG
  cerr << "\n In return list: " << endl;
#endif
           }
        )
        | ( '['
          . lparamlist
          . ']'
          )
      )
      . ([ \t]*
#          @{ buffer_.append(*p); }
        )
      :> '=' . WSOC*
    );
    # }}}2

  # a line in the function body {{{2
  funcline := |*
    ([ \t]+)
      => { cout.write(ts, te-ts); };

    ('...' . [ \t]* . EOL)
      => { cout.write(ts, te-ts); };

#    ('%' @{ tmp_p = p + 1; } . garble_comment_line);
    (comment_block)
      => {
           assert(p >= tmp_p-1);
           cout.write(tmp_p, p - tmp_p+1);
           fcall in_comment_block;
         };

    # automatically add return value fields to retval_list_
    (
     # matlab identifier (which can be a return value and a structure)
     (IDENT
        %{tmp_string.assign(ts,p-ts);})
     . '.'
     # matlab identifer (fieldname)
     . (IDENT_W_DOT >(st_tok) %{tmp_p2 = p;} )
     # if a value is assigned to this field, the field is generated/modified
     . [ \t]* . '=' . (^'=')
    )
    => {
      fhold;
      // store fieldname
      assert(tmp_p2 >= tmp_p);
      string s(tmp_p, tmp_p2 - tmp_p);
      cout << tmp_string << "." << s << "=";
      // typedef of iterators
      typedef DocuList     :: iterator list_iterator;
      typedef DocuListMap  :: iterator map_iterator;
      typedef DocuBlock    :: iterator iterator;

      // check wether first IDENT is a return value
      iterator it = find(returnlist_.begin(), returnlist_.end(), tmp_string);
      if(it != returnlist_.end())
      {
        // if it is a return value...
        // ... check wether its found field is still missing a DocuBlock in the
        // retval list.
        bool missing = true;
        map_iterator rvoit = retval_list_.find(tmp_string);
        if(rvoit != retval_list_.end())
        {
          list_iterator lit = (*rvoit).second.find(s);
          if(lit != (*rvoit).second.end())
            missing = false;
        }
        // if it is missing, add an empty docu block
        if(missing)
        {
          retval_list_[tmp_string][s] = DocuBlock();
        }
      }
    };

    # automatically add parameter fields to required_list_
    (
     # matlab identifier (which can be a parameter and a structure)
     (IDENT
        %{tmp_string.assign(ts,p-ts);})
     . '.'
     # matlab identifer (fieldname)
     . (IDENT_W_DOT
         >(st_tok)
       )
    )
    => {
      // store fieldname
      assert(p >= tmp_p);
      string s(tmp_p, p - tmp_p+1);
      cout << tmp_string << "." << s;
      typedef DocuList     :: iterator list_iterator;
      typedef DocuListMap  :: iterator map_iterator;
      typedef DocuBlock    :: iterator iterator;

      // check wether first IDENT is a parameter
      iterator it = find(paramlist_.begin(), paramlist_.end(), tmp_string);
      if(it != paramlist_.end())
      {
        // if it is a parameter ...
        // ... check wether its found field is still missing a DocuBlock in the
        // return, optional and the required list.
        bool missing = true;
        map_iterator rvoit = retval_list_.find(tmp_string);
        if(rvoit != retval_list_.end())
        {
          list_iterator lit = (*rvoit).second.find(s);
          // found match in retval list
          if(lit != (*rvoit).second.end())
            missing = false;
        }
        map_iterator moit = optional_list_.find(tmp_string);
        if(moit != optional_list_.end())
        {
          // found match in optional list
          list_iterator lit = (*moit).second.find(s);
          if(lit != (*moit).second.end())
            missing = false;
        }
        map_iterator roit = required_list_.find(tmp_string);
        if(roit != required_list_.end())
        {
          // found match in required list
          list_iterator lit = (*roit).second.find(s);
          if(lit != (*roit).second.end())
            missing = false;
        }
        // in case it IS missing, add an empty field to the required block.
        if(missing)
        {
          required_list_[tmp_string][s] = DocuBlock();
        }
      }
    };

    # add a @deprecated command to function declaration if disp_deprecated is
    # used in function body
    ('disp_deprecated' . [ \t]*
      . (
          ';'
            @{tmp_string.assign("");}
          |
          '(' . [\t ]* . "'"
          . ([^\n']*
              >(st_tok)
              %(string_tok)
            )
          . "'" . [\t ]* . ')' . [\t ]* . ';'
        )
      . [\t ]* . EOL
    )
      => {
        string s;
        if(tmp_string.empty())
        {
          s.assign("@deprecated function deprecated\n");
        }
        else
        {
          s.assign("@deprecated method deprecated, use \'" + tmp_string + "\' instead.\n");
        }
        docuextra_.push_back(s);
        fhold;
      };

    # simple matlab identifier
    (IDENT)
      => { cout.write(ts, te-ts); };

    # translate curly brackets in edgy brackets, because otherwise the doxygen
    # parser breaks.
    ('{')
      => { cout << '['; };

    ('}')
      => { cout << ']'; };

    # simply output all other characters
    (default - [\n{}])
      => { cout << fc; };

    # after EOL try to check for new function
    EOL
      => { cout << fc; fgoto funcbody; };

  *|;
  # }}}2

  # function body {{{2
  funcbody := |*

      # things that got replaced in function body {{{4
      ('% TO BE ADJUSTED TO NEW SYNTAX\n')
        => {
          new_syntax_ = true;
          cout << "*/\n"; //cout << "add to special group */\n";
        };

      # a comment block
      (comment_block)
        => {
          assert(p+1 >= tmp_p);
          cout.write(tmp_p, p - tmp_p+1);
          fcall in_comment_block;
        };

      # empty line
      ([ \t]* . EOL)
        => { cout << '\n'; };

      #}}}4

      # things that could end the function body {{{4
      # line not beginning with words 'function' or 'end'
      ([ \t]*
       . ( (default - [ \r\t\n%])+ - ('function'|'end') )
      )
        => {
          p = ts-1;
          // further parse the function body line
#ifdef DEBUG
debug_output("in funcbody: goto funcline", p);
#endif
          fgoto funcline;
        };

      # line only containing word 'end'
      # the keyword needs to be in the same indentation level as beginning function
      ([ \t]* . 'end' . ';'* . (WSOC | EOL ) )
          => {
              if(is_class_ && class_part_ == Method)
              {
                tmp_string.assign(ts,p-ts+1);

                if(tmp_string.find("e") == funcindent_)
                {
                  end_function();
#ifdef DEBUG
debug_output("in funcbody: goto methods", p);
#endif
                  fgoto methods;
                }
              }
              // else
              p=ts-1;
              // further parse the function body line
#ifdef DEBUG
debug_output("in funcbody: goto funcline 2", p);
#endif
              fgoto funcline;
          };

      # line beginning with word 'function'
      ([ \t]*. 'function ')
      {
        p = ts-1;
        // end the previous function if existent
        end_function();
#ifdef DEBUG
debug_output("in funcbody: goto main", p);
#endif
        fgoto main;
      };

      (EOF) $eof(end_of_file);

      # }}}4

  *|;
   # }}}2

  # fill a docublock list with input {{{2
  fill_list := |*

  # match an argument
  ( doc_begin . [ \t]*
    . "'"? . ( ([A-Za-z][A-Za-z0-9_{},()[\].]*)  >{tmp_p3 = p;} %{tmp_p2 = p;} ) . "'"? . [ \t]* . ":" @(st_tok)
    . ( default - '\n' )* . EOL
  )
    => {
      assert(tmp_p2 >= tmp_p3);
      tmp_string.assign(tmp_p3, tmp_p2 - tmp_p3);
      //    std::cout << tmp_string << '\n';
      assert(p >= tmp_p);
      (*clist_)[tmp_string].push_back(string(tmp_p+1, p - tmp_p));
    };

  # expand the paragraph for last argument matched
  ( doc_begin . [ \t]*
    # at least one word (non white-space characters and no double-colon)
    . ( default - [ \r\t:\n] )+ .
    # followed by something that is a white-space or a new-line, i.e *no*
    # double-colon
    (
     EOL
     |
     [ \t]+ . (EOL | [^ \r\n\t:] . (default - '\n')* . EOL)
     # [ \t] . (default - '\n')* . EOL
    )
  )
    => {
      assert(p+1 >= tmp_p);
      string s(tmp_p, p - tmp_p + 1);
      (*clist_)[tmp_string].push_back(s);
      /*cout << "add something results in\n" << (*clist_)[tmp_string];*/
    };

  # return on empty line
  ( doc_begin . [ \t]* . EOL )
    => { /*cout << "empty line\n";*/ fret; };

   # end of comment block
  ( [\t ]* . ( (default - '%') | EOL) )
    => {
      p =ts-1;
      // cout << "*/\n";
      fret;
    };

  *|; #}}}2

  # parse body of documentation block {{{2
  doxy_get_body := |*

    # special lists {{{4

    # begin required_list
    ( doc_begin . [ \t]*
      . /required fields of /i
      . (IDENT >(st_tok) %(string_tok) )
      . [ \t]* . ':' . [ \t]* . EOL
    )
      => {
        //cout << tmp_string << '\n';
        clist_ = &(required_list_[tmp_string]);
        docline = false;
        fcall fill_list;
      };

    # begin optional_list
    ( doc_begin . [ \t]*
      . /optional fields of /i
      . (IDENT
          >(st_tok)
          %(string_tok) )
      . [ \t]* . ':' . [ \t]* . EOL )
      => {
        clist_ = &(optional_list_[tmp_string]);
        docline = false;
        fcall fill_list;
      };

    # begin optional_list
    ( doc_begin . [ \t]*
      . /generated fields of /i
      . (IDENT
          >(st_tok)
          %(string_tok) )
      . [ \t]* . ':' . [ \t]* . EOL )
      => {
        clist_ = &(retval_list_[tmp_string]);
        docline = false;
        fcall fill_list;
      };

    # begin parameter list
    ( doc_begin . [ \t]*
      . /parameters/i . [ \t]* . ':'
      . [ \t]* . EOL )
      => {
        clist_ = &param_list_;
        docline = false;
        fcall fill_list;
      };

    # begin return list
    ( doc_begin . [ \t]*
      . /return values/i . [ \t]* . ':'
      . [ \t]* . EOL )
      => {
        clist_ = &return_list_;
        docline = false;
        fcall fill_list;
      };
    #}}}4

    # default substitutions {{{4

    # empty line
    ( doc_begin . [ \t]* . EOL )
      => {
        /*cout << "*\n  ";*/
        docubody_.push_back("\n");
        docline = false;
      };

    # paragraph line
    ( [ \t]* . '%' )
      => {
        if(!docline)
        {
          docline = true;
          tmp_p = p;
        }
      };

    # paragraph line with "see also" substituted by "@sa"
    ( /see also/i . ':'? )
      => {
        string s;
        assert(ts > tmp_p);
        s.assign(tmp_p+1, ts - tmp_p-1);
        docubody_.push_back(s+"@sa");
        tmp_p = p;
      };

    # lines that could end doxyblock {{{6
    # words
    #    ( default - [ \t:%'`\n] )+
    ( default - [ \t:%\r\n] )+ @(end_doxy_block);

    # non-words/non-whitespace
    #    ([:'`]) => {
    (':') @(end_doxy_block) ;


    # whitespace only
    ( [ \t] );

    # titled paragraph
    ( ':' . EOL )
      @(end_doxy_block)
      @{ if(docline)
         {
           assert(ts > tmp_p);
           docubody_.push_back("@par " + string(tmp_p+1, ts - tmp_p-1)+"\n");
           docline = false;
         }
       };
    # }}}6
    # }}}4

    # end of line {{{4
    ( EOL )
       @(end_doxy_block)
       @{ if(docline)
          {
            int offset = ( latex_begin ? 0 : 1 );
            assert(p >= tmp_p + offset);
            docubody_.push_back(string(tmp_p+1, p - tmp_p - offset));
            docline = false;
          }
        };
      # }}}4

  *|;
  #}}}2

  # doxy header parsing {{{2
  # swallow the synopsis line
  doxyfunction_garble := |*
    garbage = ( (default - '\n' )* -- '...' );

    ( doc_begin . (garbage . '...')+ . [\t ]* .  EOL );

    ( doc_begin . (garbage . '...')* . garbage . EOL )
      => { fgoto doxy_get_brief; };
  *|;


  # read first paragraph
  doxy_get_brief := |*

    # read in one comment line
    ( doc_begin . [\t ]*
      . (default - [\r\n\t ]) . (default - '\n')* . EOL
    )
      => {
        /* cout << "*"; cout.write(tmp_p, p - tmp_p+1); */
        assert(p >= tmp_p);
        docuheader_.push_back(string(tmp_p, p - tmp_p+1));
      };

    # empty line
    ( doc_begin . [\t ]* . EOL )
      => {
        /*cout << "*\n";*/
#ifdef DEBUG
  debug_output("in doxy_get_brief: goto: doxy_get_body", p);
#endif
        fgoto doxy_get_body;
      };

    # end of comment block;
    ( [\t ]* . [^%] )
      => {
        p=ts-1;
#ifdef DEBUG
   debug_output("in doxy_get_brief: end!!", p);
#endif
        //cout << "*/\n";
        if(is_class_)
        {
#ifdef DEBUG
  debug_output(" in_doxy_get_brief: this is a class",p);
#endif

          if(class_part_ == Header)
          {
#ifdef DEBUG
  debug_output("  in_doxy_get_brief: method: goto classbody",p);
#endif
            end_of_class_doc();
            fgoto classbody;
          } else if(class_part_ == Method || class_part_ == AtMethod)
          {
#ifdef DEBUG
  debug_output("  in_doxy_get_brief: method: goto funcbody",p);
#endif
            print_function_synopsis();
            fgoto funcbody;
          }
          else if(class_part_ == MethodDeclaration)
          {
#ifdef DEBUG
  debug_output("  in_doxy_get_brief: method: goto funcdef",p);
#endif
            fgoto funcdef;
          }
          else if(class_part_ == Property)
          {
#ifdef DEBUG
  debug_output("  in_doxy_get_brief: method: goto propertybody",p);
#endif
            fgoto propertybody;
          }
          else if(class_part_ == InClassComment)
          {
            class_part_ = Method;
            fgoto methods;
          }
        }
        else
        {
          print_function_synopsis();
          fgoto funcbody;
        }
      };

  *|;
  # }}}2

  # garble synopsis line and then parse the documentation header {{{2
  doxyheader := (
    '%' . [ \t]* .
       (
        ('function '|'classdef ') @{ p = tmp_p-2; fgoto doxyfunction_garble; }
       )
      $!{
#ifdef DEBUG
        debug_output("doxy_get_brief",p);
#endif
        p = tmp_p - 2;
        fgoto doxy_get_brief;
      }
   ); #}}}2

  # helper for setting the access specifier {{{2
  paramaccess =
    ( ('SetAccess' . WSOC* . '=' . WSOC*
      . ( (/public/i
            @{ access_.full = Public;
               access_.set = Public;
             } )
        | ( /protected/i
            @{ access_.full =
                 (access_.get == Public ? Public : Protected );
               access_.set = Protected;
             } )
        | ( /private/i
            @{ access_.full = access_.get;
               access_.set = Private;
             } )
        )
      )
     | ( 'GetAccess' . WSOC* . '=' . WSOC*
      . ( ( /public/i
            @{ access_.full = Public;
               access_.get = Public;
             } )
        | ( /protected/i
            @{ access_.full =
                 (access_.set == Public ? Public : Protected );
               access_.get = Protected;
             } )
        | ( /private/i
            @{ access_.full = access_.set;
               access_.get = Private;
             } )
        )
       )
     | ( 'Access' . WSOC* . '=' . WSOC*
      . ( ( /public/i
            @{ access_.full = Public;
               access_.get = Public;
               access_.set = Public;
             } )
        | ( /protected/i
            @{ access_.full = Protected;
               access_.get = Protected;
               access_.set = Protected;
             } )
        | ( /private/i
            @{ access_.full = Private;
               access_.get = Private;
               access_.set = Private;
             } )
        )
       )
      ); #}}}2

  # method and property params {{{2
  methodparam =
   (
    ( paramaccess )
    | ( ( 'Abstract' . [^,)]* )
        @{
           methodparams_.abstr = true;
         } )
    | ( ( 'Static' . [^,)]* )
        @{
           methodparams_.statical = true;
         } )
    | ( ('Hidden'|'Sealed')
         . [^,)]* )
   );

  propertyparam =
   (
    ( paramaccess )
    | ( ( 'Constant' . [^,)]* )
        @{
           propertyparams_.constant = true;
         } )
    | ( ('Transient'|'Dependent'|'Hidden'|'SetObservable'|'Abstract')
         . [^,)]* )
   );

  methodparams =
   (
    '(' . WSOC*
    . methodparam
    . ( WSOC* . ',' . WSOC* . methodparam )* . WSOC* . ')'
   );

  propertyparams =
   (
    '(' . WSOC*
    . propertyparam
    . ( WSOC* . ',' . WSOC* . propertyparam )* . WSOC* . ')'
   ); #}}}2

  # swallowing events {{{2
  events := (
    ( ([ \t]* .  'e') >st_tok
      . 'nd' . [ \t;]* . EOL
        @{ tmp_string.assign(tmp_p, p - tmp_p);
                if(tmp_string.find("e") == eventindent_)
                  {
                  fgoto classbody;
                  }
              }
      | (default* - 'end') . EOL )*
   ); #}}}2

  # methods and properties {{{2
  # methods {{{4
  methods := |*
# kommentare, newlines
# nur bei keyword 'function' => goto funcdef
# abstrakter fall, eine weitere Regel wird benötigt.
# end => classbody
    (empty_line) => {
      end_method();
      cout << "\n";
    };

    ([ \t]* . 'function' )
      => {
        tmp_string.assign(ts, te - ts+1);
        funcindent_ = tmp_string.find_first_not_of(" \t");
#if DEBUG
    {
      ostringstream oss;
      oss << "in methods: funcindent: " << funcindent_;
      debug_output(oss.str(), p);
    }
#endif
        p=ts+funcindent_-1;
        fgoto funct;
       };

    ([ \t]* . 'end' . [ \t;]* . ('%' . garble_comment_line_wo_eol)? . EOL )
      => {
           end_method();
#if DEBUG
    debug_output("in methods: found end keyword, goto classbody",p);
#endif
           fgoto classbody;
         };

    ([ \t]* . '%' ) => {
#if DEBUG
    debug_output("in methods: garble comment line",p);
#endif

      p = ts-1;
      class_part_ = InClassComment;
/*      fcall in_comment_block; */
      fgoto expect_doxyblock;
    };

    # if we reach this: method declaration without definition is found
    ([ \t]* . [^% \t\n]) =>
    {
#if DEBUG
    debug_output("in methods: found method declaration, going to funcdef",p);
#endif
      class_part_ = MethodDeclaration;
      p = ts-1;
      fgoto funcdef;
    };

      *|;


  methodsheader := (
    [ \t]* . methodparams? . [ \t;]* . ( '%' . garble_comment_line_wo_eol )? . EOL
         @{
            print_access_specifier(access_.full);
            fgoto methods;
          }
              );

   #}}}4

  # single property {{{4
  prop = ( ( [ \t]* . (IDENT) >st_tok
                 %{
            string s(tmp_p, p - tmp_p);
            property_list_.push_back(s);
//            cout << propertyparams_.ccprefix() << " " << s;
            }
          )
        . ( ';' @{defaultprop_ = "";}
            |
            ([ =]+ %st_tok . ('[' . [^\]]* . ']')? . [^;[]* .';')
            @{
              defaultprop_ = string(tmp_p, p - tmp_p);
             }
          )
            . [ \t]* .
        ( '%' %st_tok . ( default - [\r\n] )*
          . EOL
            @{
               docuheader_.push_back(string(tmp_p, p - tmp_p+1));
               end_of_property_doc();
             } | EOL @{ end_of_property_doc(); }
        )
      );

  #}}}4

  property := ( prop* );

  #property body {{{4
  propertybody = (
    ( [ \t]* . ( 'end' . [ \t;]* ) . ('%' . garble_comment_line_wo_eol )? . EOL )
      @{
         fgoto classbody;
       }
    |
    (prop)
    |
    ( (empty_line) @{ cout << "\n";} )
    |
    ( ([ \t]* . '%') @{ fhold; fgoto expect_doxyblock; } )
    );

  properties := ( (
    WSOC* . propertyparams? . [ \t;]* . ('%' . garble_comment_line_wo_eol )? . EOL @{
        print_access_specifier(access_.full);
        }
    . propertybody* )
      );
  #}}}4
  #}}}2

  # class body {{{2
  classbody := |*

    # a comment block
    (comment_block)
      => {
        cout.write(tmp_p, p - tmp_p+1);
        fcall in_comment_block;
      };

    (WSOC) => { cout.write(ts, te-ts); };

    (EOL) => { cout << "\n"; };

    ('end' . [ \t]* ';'?) => {
      cout << "\n};\n";
      for(  list<string>::iterator it = namespaces_.begin();
            it != namespaces_.end(); ++it)
      {
        cout << "}\n";
      }
    };

    ([ \t]* . 'properties')
      => {
        propertyparams_ = PropParams();
        access_ = AccessStruct();
        class_part_ = Property;
        fgoto properties;
      };
    ([ \t]* . 'methods')
      => {
        methodparams_ = MethodParams();
        access_ = AccessStruct();
        class_part_ = Method;
        fgoto methodsheader;
      };
    ([ \t]* . 'events')
      => {
        std::string tmp_string(ts, te-ts);
        eventindent_ = tmp_string.find("e");
        class_part_ = Event;
        fgoto events;
      };
  *|; #}}}2

  # doxyblock expect {{{2
  # after function declaration expect a documentation block or the function
  # body
  expect_doxyblock :=
  (
    doc_begin
      @{
        //cout << "/*";
        p--;
        fgoto doxyheader;
      }
  )
 $!{
    fhold;
#ifdef DEBUG
    debug_output("stopping expect_doxyblock", p);
#endif
    if(is_class_)
    {
      if(class_part_ == Header)
      {
        end_of_class_doc();
        fgoto classbody;
      } else if(class_part_ == Method || class_part_ == AtMethod)
      {
        string endstringtest;
        endstringtest.assign(p, 100);
        string::size_type first_char = endstringtest.find_first_not_of(" \t");
        if (endstringtest.substr(first_char, 3) == "end")
        {
          p += first_char+4;
          print_function_synopsis();
          end_function();
          fgoto methods;
        }
        else
        {
          print_function_synopsis();
          fgoto funcbody;
        }
      }
      else if(class_part_ == Property)
      {
        fgoto propertybody;
      }
      else if(class_part_ == InClassComment || class_part_ == MethodDeclaration)
      {
        class_part_ = Method;
        fgoto methods;
      }
      else{
        cerr << "Do not know where to go from here. Classpart " << ClassPartNames[class_part_] << " is not handled.\n";
      }
    }
    else
    {
      print_function_synopsis();
      fgoto funcbody;
    }
  };
  #}}}2

  # function declaration {{{2
  funcdef = (
      (WSOC)* .
      # return values (if found opt = true)
      (lparams)? .
      # matlab identifier (function name stored in cfuncname_)
      ( ('get.' @{is_getter_ = true;} | 'set.' @{is_setter_=true;} )? .
        IDENT
          >st_tok
          %{
            cfuncname_.assign(tmp_p, p - tmp_p);
#ifdef DEBUG
  cerr << "\n Identifier of function: " << cfuncname_ << endl;
#endif
            is_script_ = false;
          }
      )
      . WSOC*
      . (
           '('
           # parameter list
           . ( paramlist
               %{
                 if(paramlist_.size() == 1 && paramlist_[0] == "this")
                 { paramlist_.clear(); }
               }
             )
           . ')' . ( [ \t] | ('%' @{ tmp_p=p; comment_found=true; }
             . garble_comment_line_wo_eol) | ( ';' ) )*
           . EOL
           @{
             if(comment_found)
             {
                tmp_string.assign(tmp_p+1, p - tmp_p-1);
                tmp_string = string("/* ") + tmp_string + string("*/");
             }
             else
             {
                tmp_string = "";
             }
             comment_found = false;
             if(is_class_ && class_part_ == MethodDeclaration )
             {
               class_part_ = Method;
#if DEBUG
    debug_output("in funcdef: end of method declaration, returning to methods",p);
#endif
               fgoto methods;
             }
             else
             {
//               cout << tmp_string << "{\n";
               // check for documentation block
               fgoto expect_doxyblock;
             }
           }
        # no parameter list && first function => script
        | (( [ \t]
              |
             ('%' @{ tmp_p=p; comment_found=true; }
              . garble_comment_line_wo_eol) | ( ';' ) )* . EOL)
            @{
                 if(comment_found)
                 {
                   tmp_string.assign(tmp_p+1, p - tmp_p-1);
                   tmp_string = string("/* ") + tmp_string + string("*/");
                 }
                 else
                 {
                   tmp_string = "";
                 }
                 comment_found = false;
#if DEBUG
  debug_output("in funcdef: script && no parameters: expect doxyblock",p);
#endif
                 if(is_class_ && class_part_ == MethodDeclaration)
                 {
                    class_part_ = Method;
                    fgoto methods;
                 }
                 else
                 {
                   fgoto expect_doxyblock;
                 }
             }
        )
      );

  funct :=
  (
    (
      comment_block @in_c_block
      | [ \t]*. EOL
    )*
    . [\t]* . 'function'
    . funcdef
  ) $eof( end_of_file ) ; #}}}2

  # no function definition => a script {{{2
  script := (default)
    @{
       string :: size_type found = filename_.rfind("/");
       if(found == string :: npos)
         found = -1;
       string funcname = filename_.substr(found+1, filename_.size()-3-found);
       cfuncname_.assign( funcname );
  /*     cout << "noret::substitute ";
       if(!is_first_function_)
         cout << "mtoc_subst_" << fnname_ << "_tsbus_cotm_";
       cout << funcname << "() {\n";*/
       is_script_ = true;
       fhold;
       fgoto expect_doxyblock;
     }; #}}}2

  # class definitions {{{2
  classparams =
      '(' . [^)]* . ')';

  superclass =
    ( ( IDENT_W_DOT ) >{ cout << "public ::"; }
                      @{ if(*p == '.')
                           cout << "::";
                         else cout << *p; } );

  superclasses = (
      '<' @{ cout << "\n  :"; } . WSOC* . superclass . WSOC*
          . ('&' @{ cout << ",\n   "; } . WSOC* . superclass . WSOC*)* );

  classdef = (
      'classdef' . WSOC* .
      # matlab identifier (class name stored in classname_)
      ( IDENT
          >st_tok
          %{
            classname_.assign(tmp_p, p - tmp_p);
            is_class_ = true;
            cout << "class " << classname_;
          }
      )
      . WSOC*
      . classparams?
      . WSOC*
      . superclasses?
      . [ \t;]*
      . ( '%'. garble_comment_line_wo_eol )?
      EOL
      @{
        cout << " {\n";
        fgoto expect_doxyblock;
      } );

  class :=
  (
    (
      comment_block @in_c_block
      | [ \t]*. EOL
    )*
    . [\t]* . 'classdef'
    . classdef
  ) ;
  # }}}2

  # main loop {{{2
  expect_function_script_or_class =
  (
    # either we find a function or classdef definition with a possibly
    # preceding comment block or we have a script
    ( any @{ fhold; tmp_p = p; } .
    (
      [ \t]*. '%' . (any - '\n')* . EOL
      | [ \t]*. EOL
    )*
    . [\t]*
    . ( 'function' @{
                     p=tmp_p;
                     if(is_class_ && class_part_ == Header)
                       class_part_ = AtMethod;
                     fgoto funct;
                    }
      | 'classdef' @{
                     p=tmp_p;
                     fgoto classdef;
                    }
      ) )
  $!{
#ifdef DEBUG
    debug_output("goto script",p);
#endif
    p=tmp_p;
    fgoto script;
  }
  );

  main := expect_function_script_or_class*;
  # }}}2

}%%

void MFileScanner :: print_pure_function_synopsis()
{
  // do we have a constructor?
  if(is_class_ && (cfuncname_ == classname_))
    returnlist_.clear();
  else{
    if(returnlist_.size() == 0)
      cout << "noret::substitute ";
    else if(returnlist_.size() == 1)
      cout << "ret::substitutestart::" << returnlist_[0] << "::retsubstituteend ";
    else
    {
      cout << "rets::substitutestart::";
      for(unsigned int i=0; i < returnlist_.size(); ++i)
      {
        cout << returnlist_[i] << "::";
      }
      cout << "retssubstituteend ";
    }
  }

  bool first = true;
  if(is_first_function_)
  {
    if(is_class_ && class_part_ == AtMethod)
      cout << namespace_string() << classname_ << "::";
  }
  else
    cout << "mtoc_subst_" << fnname_ << "_tsbus_cotm_";

  cout << cfuncname_;

  if(paramlist_.size() == 0)
    cout << "()\n  ";
  else
  {
#if DEBUG
    cerr << "paramlist size of " << cfuncname_ << ": " << paramlist_.size() << " first element: " << paramlist_[0] << endl;
#endif
    cout << "(";
    for(unsigned int i=0; i < paramlist_.size(); ++i)
    {
      if(!first)
        cout << ",";
      else
        first = false;

      std::string typen;// = "matlabtypesubstitute";
      get_typename(paramlist_[i], typen);
      cout << typen << " " << paramlist_[i];
    }
    cout << ")";
  }
}

void MFileScanner :: print_function_synopsis()
{
  if(is_getter_ || is_setter_)
  {
   cout << "/* \n";
  }
  if(is_class_ && (class_part_ == Method || class_part_ == MethodDeclaration) )
    cout << methodparams_.ccprefix();

  print_pure_function_synopsis();

  if(is_class_ && class_part_ == MethodDeclaration )
    cout << methodparams_.ccpostfix() << "\n";
  else
    cout << " {\n";
}

void MFileScanner :: print_access_specifier(AccessEnum & access)
{
  if(access == Public)
    cout << "public:\n";
  else if(access == Protected)
    cout << "protected:\n";
  else if(access == Private)
    cout << "private:\n";
}

// constructor
MFileScanner :: MFileScanner(istream & fin, const std::string & filename,
                             const std::string & conffilename, bool latex_output) :
  fin_(fin), filename_(filename),
  latex_output_(latex_output), cscan_(filename_, conffilename),
  fnname_(filename), namespaces_(),
  line(1),
  ts(0), have(0), top(0),
  opt(false), new_syntax_(false),
  is_script_(false), is_first_function_(true),
  is_class_(false), is_setter_(false), is_getter_(false),
  classname_(), funcindent_(0), eventindent_(0),
  class_part_(Header),
  access_(), propertyparams_(), methodparams_(), property_list_()
{
  string::size_type found = fnname_.find_last_of('/');
  std::string dirname;
  if(found != string::npos)
    dirname = filename.substr(0, found);

  list<string> namespaces;
  string classname;
  string::size_type enddir = dirname.size();
  string::size_type ppos = 0;
  while (ppos != string::npos)
  {
    ppos = dirname.find_last_of('/', enddir);
    string directory;
    if(ppos == string::npos)
      directory = dirname.substr(0, enddir+1);
    else
      directory = dirname.substr(ppos+1, enddir-ppos);

    if(directory[0] == '+')
    {
      namespaces_.push_front(directory.substr(1));
    }
    else if(directory[0] == '@')
    {
      classname_ = directory.substr(1);
      is_class_ = true;
      if(classname_
          != fnname_.substr(fnname_.find_last_of('/')+1, classname_.size()))
      {
        class_part_ = AtMethod;
        cout << "#include \"" << classname_ << ".m\"" << endl;
      }
    }
    else
      break;
    enddir = ppos - 1;
  }
  for (list<string>::iterator it = namespaces_.begin();
       it != namespaces_.end(); ++it)
    cout << "namespace " << *it << "{\n";

  found = fnname_.rfind("/");
  if(found != string::npos)
    fnname_ = fnname_.substr(found+1);
  for( std::string::size_type i = 0; i < fnname_.size(); ++i )
  {
    if(fnname_[i] == '@')
      fnname_[i] = '_';
    else if(fnname_[i] == '.')
      fnname_[i] = '_';
  }

  cscan_.execute();
};

// run the scanner
int MFileScanner :: execute()
{
  std::ios::sync_with_stdio(false);

  %% write init;

  /* Do the first read. */
  bool done = false;
  while ( !done )
  {
    char *p = buf + have;
    char *tmp_p = p, *tmp_p2 = p, *tmp_p3 = p;
    string tmp_string;
    bool docline = false;
    bool latex_begin = true;
    bool comment_found = false;
    int space = BUFSIZE - have;

    if ( space == 0 )
    {
      /* We filled up the buffer trying to scan a token. */
      cerr << "OUT OF BUFFER SPACE" << endl;
      exit(-1);
    }

    fin_.read( p, space );
    int len = fin_.gcount();
    char *pe = p + len;
    char *rpe = pe;
    char *eof = 0;

    /* If we see eof then append the EOF char. */
    if ( fin_.eof() )
    {
      char eof_c = *pe;
      *pe = '\n';
      pe++;
      *pe = eof_c;
      eof = pe;
      rpe = pe;

      done = true;
    }
    else
    {
      /* Find the last newline by searching backwards. This is where
       * we will stop processing on this iteration. */
      while ( pe >= p )
      {
        if( *pe != '\n')
          pe--;
        else
        {
          if(pe >= p+3
              && *(pe-1) == '.' && *(pe-2) == '.' && *(pe-3) == '.')
            pe-=3;
          else
            break;
        }
      }
    }

    %% write exec;

    /* Check if we failed. */
    if ( cs == MFileScanner_error )
    {
      /* Machine failed before finding a token. */
      cerr << std::string(filename_) << ": PARSE ERROR in line " << line << endl;
      debug_output("Grrrr!!!!", p);
      exit(-1);
    }

    /* Now set up the prefix. */
    if ( ts == 0 )
    {
      have = rpe - pe;
      /* cerr << "memmove by " << have << "bytes\n";*/
      memmove( buf, pe, have );
    }
    else
    {
      have = rpe - ts;
      /* cerr << "memmove by " << have << "bytes to ts\n";*/
      memmove( buf, ts, have );
    }

    if ( ts != 0 )
    {
      te -= (ts-buf);
      ts = buf;
    }
  }

  return 0;
}

// escape '@' and '\' characters in string \a s
const string & MFileScanner::escape_chars(std::string & s)
{
  string::size_type found = s.find_first_of("@\\");
  while(found != string::npos )
  {
    s.insert(found, "\\");
    found = s.find_first_of("@\\",found+2);
  }
  return s;
}

// standard brief text (replace '_' -> ' ' in s)
const string & MFileScanner::replace_underscore(std::string & s)
{
  string::size_type found = s.find("_");
  while(found != string::npos )
  {
    s[found] = ' ';
    found = s.find("_", found+1);
  }
  return s;
}

// pretty print the documentation block \a block
void MFileScanner::write_docu_block(const DocuBlock & block)
{
  bool add_prefix   = false;
  bool latex_begin  = true;
  bool not_verbatim = true;
  for( unsigned int i = 0; i < block.size(); i += 1 )
  {
    // begin all documentation lines after the first one with an asterisk
    if(add_prefix)
      cout << "* ";

    add_prefix = false;
    // read in new line of docu block
    const string & s = block[i];

    // parse for special comments
    string::size_type j=0;
    const char * tokens = "\'`@\n";
    bool last_char_escaped = false;
    for( string::size_type i = 0; j < s.size(); i=j )
    {
      j=s.find_first_of(tokens,i+1);
      if(j==string::npos)
        j=s.size();
      if(s[j-1] == '\\' && not_verbatim && latex_begin)
        --j;
      // respect @code and @verbatim blocks
      if(s[i] == '@')
      {
        if(s.substr(i+1,4) == "code" || s.substr(i+1,8) == "verbatim")
          not_verbatim = false;
        else if(s.substr(i+1,7) == "endcode" || s.substr(i+1,11) == "endverbatim")
          not_verbatim = true;
        cout << s.substr(i,j-i);
      }
      // use typewriter fonts for words in single quotes
      else if(s[i] == '\'' && not_verbatim && latex_begin)
      {
        if(j != s.size() && s[j] == '\'' && !last_char_escaped)
        {
          cout << "<tt>" << s.substr(i+1, j-i-1) << "</tt>";
          ++j;
        }
        else
          cout << s.substr(i,j-i);
      }
      // use latex output for words in backtick quotes
      else if(s[i] == '`' && not_verbatim)
      {
        string lout;
        if(!last_char_escaped)
        {
          // in case of double backtick quotes, use latex block
          if(s[i+1] == '`')
          {
            if(latex_begin)
              lout = "@f[";
            else
              lout = "@f]";
            ++i;
            j=s.find_first_of(tokens,i+1);
            if(j==string::npos)
              j=s.size();
          }
          else
            lout = "@f$";
          if(latex_begin)
            latex_begin = false;
          else
            latex_begin = true;
          ++i;
        }
        else
        {
          lout = "";
        }
        cout << lout << s.substr(i, j-i);
      }
      // new line
      else if(s[i] == '\n')
      {
        cout << "\n  ";
        if(latex_begin)
          add_prefix = true;
        else
        {
          cout << "  ";
          add_prefix = false;
        }
      }
      else
      {
        cout << s.substr(i,j-i);
      }
      if(s[j-1] != '\\' && s[j] == '\\')
      {
        last_char_escaped = true;
      }
      else
        last_char_escaped = false;
      if(s[j] == '\\')
        ++j;
    }
  }
}

// pretty print the documentation block list \a list for the list item named \a
// item_text. If docu blocks are empty, \a alternative is used. The alternative
// is normally read in by the confscanner.
void MFileScanner::write_docu_list(const DocuList & list,
                                   const string & item_text,
                                   const DocuList & alternative,
                                   const string separator = string())
{
  typedef DocuList :: const_iterator list_iterator;
  list_iterator lit = list.begin();
  // iterate over documentation blocks
  for(; lit != list.end(); ++lit)
  {
    cout << "* " << item_text << " " << (*lit).first << separator << "    ";
    const DocuBlock & block = (*lit).second;
    // block is empty?
    if(block.empty())
    {
      // then look for alternative documentation block form global
      // configuration file or use default text generated from variable name.
      list_iterator alit = alternative.find((*lit).first);
      if(alit == alternative.end() || (*alit).second.empty())
      {
        std::string s((*lit).first);
        cout << replace_underscore(s) << "\n  ";
      }
      else
        write_docu_block((*alit).second);
    }
    else
      write_docu_block(block);
  }
}

// pretty print a documentation block list map \a listmap with prepended title
// \a text. If listmap entry is empty, \a altlistmap is used instead.
void MFileScanner::write_docu_listmap(const DocuListMap & listmap,
                                      const string & text,
                                      const DocuListMap & altlistmap)
{
  typedef DocuListMap :: const_iterator map_iterator;
  if(!listmap.empty())
  {
    map_iterator mit = listmap.begin();
    for(; mit != listmap.end(); ++mit)
    {
      cout << "*\n  ";
      cout << "* " << text << (*mit).first << ":\n  ";
      map_iterator amit = altlistmap.find((*mit).first);
      write_docu_list((*mit).second, "@arg \\c", ( amit != altlistmap.end() ? (*amit).second : DocuList() ), "&nbsp;&mdash;&nbsp;");
    }
//    cout << "* </TABLE>\n  ";

  }
}

string MFileScanner::namespace_string()
{
  ostringstream oss;
  oss << "";
  for( list<string>::iterator it = namespaces_.begin();
       it != namespaces_.end(); ++it)
  {
    oss << *it << "::";
  }
  return oss.str();
}

void MFileScanner::end_of_class_doc()
{
  cout << "/** @class \"" << namespace_string() << classname_ << "\"\n  ";

  cout_ingroup();

  cout << "* @brief ";
  cout_docuheader(cfuncname_);
  cout << "*\n  ";
  cout_docubody();
  cout << "*\n ";
  cout_docuextra();
  cout << "*/\n";
}

void MFileScanner::end_of_property_doc()
{
  typedef DocuBlock :: iterator                                      DBIt;
  string typen;
  for(DBIt dit = docuheader_.begin(); dit != docuheader_.end(); ++dit)
  {
    std::string line = *dit;
    size_t found = line.find("of type");
    if(found != std::string::npos)
    {
      size_t typenstart=found+1+string("of type").length();
      size_t typenend =
        line.find_first_of( " \0", typenstart );
      typen = line.substr(typenstart, typenend - typenstart);
//      line.erase(found, typenend - found);
      line.insert(typenend, string(" \"") + typen + "\" @endlink");
      line.insert(typenstart, string("@link "));

      *dit = line;
    }
  }

  if(typen.empty())
    typen = "matlabtypesubstitute";

  cout << propertyparams_.ccprefix() << typen << " " << property_list_.back();
  if(defaultprop_.empty())
    cout << ";\n";
  else
    cout << " = " << defaultprop_ << ";\n";

  cout << "/** @var " << property_list_.back() << "\n  ";
  cout << "* @brief ";
  cout_docuheader(property_list_.back());
  cout << "*\n  ";
  cout_docubody();
  cout << "*\n ";
  cout_docuextra();
  cout << "*/\n";
  docuheader_.clear();
  docubody_.clear();
  docuextra_.clear();
}

void MFileScanner::cout_docuheader(string altheader)
{
  if(docuheader_.empty() && cscan_.docuheader_.empty())
  {
    cout << replace_underscore(altheader) << "\n  ";
  }
  else
  {
    if(! docuheader_.empty())
    {
      write_docu_block(docuheader_);
    }
    if(! cscan_.docuheader_.empty())
    {
      write_docu_block(cscan_.docuheader_);
    }
  }
  docuheader_.clear();
}

void MFileScanner :: cout_docubody()
{
  if(!docubody_.empty())
  {
    cout << "*\n  * ";
    write_docu_block(docubody_);
  }
  docubody_.clear();
  if(!cscan_.docubody_.empty())
  {
    cout << "*\n  * ";
    write_docu_block(cscan_.docubody_);
  }
}

void MFileScanner :: cout_docuextra()
{
  if(! cscan_.docuextra_.empty())
  {
    cout << "*\n  * ";
    write_docu_block(cscan_.docuextra_);
  }
  docuextra_.clear();
}

void MFileScanner :: cout_ingroup()
{
  typedef GroupSet     :: iterator group_iterator;
  // add @ingroup commands from the configuration file
  if((! groupset_.empty() || ! cscan_.groupset_.empty() ))
  {
    cout << "* @ingroup ";
    bool not_first = false;
    group_iterator git = cscan_.groupset_.begin();
    for(; git != cscan_.groupset_.end(); ++git)
    {
      if(not_first)
        cout << " ";
      else
        not_first = true;

      cout << *git;
    }
    groupset_.clear();
    cout << "\n  ";
  }
}

void MFileScanner::clear_lists()
{
#ifdef DEBUG
  std::cerr << "clear lists" << endl;
#endif
  paramlist_.clear();
  returnlist_.clear();
  param_list_.clear();
  return_list_.clear();
  required_list_.clear();
  optional_list_.clear();
  retval_list_.clear();
}

/* we come here, from an empty line in a methods block or the end of a
 * methods block
 */
void MFileScanner::end_method()
{
  if (!cfuncname_.empty())
  {
    class_part_ = MethodDeclaration;
    print_function_synopsis();
    class_part_ = Method;

    // for abstract methods: print out documentation of the abstract method
    // declaration
    if(methodparams_.abstr)
      end_function();
    else
    // otherwise: all the following comments are not related to this function
    // anymore, so we delete traces of the method name...
    {
      cfuncname_.clear();
      clear_lists();
    }
  }
  // free documentation block variables
  docuheader_.clear();
  docubody_.clear();
  docuextra_.clear();
}

void MFileScanner::get_typename(const std::string & paramname, std::string & typen)
{
  typedef DocuList :: iterator                                       DLIt;
  typedef DocuBlock :: iterator                                      DBIt;
  DLIt it  = param_list_.find(paramname);
  DocuBlock * pdb;
  if(it != param_list_.end() && !(it->second).empty())
    pdb   = &(it->second);
  else
  {
    it = return_list_.find(paramname);
    if(it != return_list_.end() && !(it->second).empty())
      pdb   = &(it->second);
    else
    {
      it = cscan_.param_list_.find(paramname);
      if(it != cscan_.param_list_.end() && !(it->second).empty())
        pdb   = &(it->second);
      else
      {
        it = cscan_.return_list_.find(paramname);
        if(it != cscan_.return_list_.end() && !(it->second).empty())
          pdb   = &(it->second);
        else
        {
          typen="matlabtypesubstitute";
          return;
        }
      }
    }
  }

  DocuBlock & db = *pdb;
  for(DBIt dit = db.begin(); dit != db.end(); ++dit)
  {
    std::string line = *dit;
    size_t found = line.find("of type");
    if(found != std::string::npos)
    {
      size_t typenstart=found+1+string("of type").length();
      size_t typenend =
        line.find_first_of( " \0", typenstart );
      typen = line.substr(typenstart, typenend - typenstart);
      line.erase(found, typenend - found);
    }
  }

  if(typen.empty())
    typen = "matlabtypesubstitute";
}

// end a function and pretty print the documentation for this function
void MFileScanner::end_function()
{
  bool is_constructor = false;
  bool is_method = false;
  if(is_class_)
  {
    if(class_part_ == Property)
      return;
    if(cfuncname_ == classname_)
      is_constructor = true;
    if(class_part_ == Method)
      is_method = true;
  }
  // end function
  if(!is_method || !methodparams_.abstr)
    cout << "}\n";
  if(is_getter_ || is_setter_)
    cout << "*/\n";
  // is the first function?
  if(is_first_function_)
  {
    if(! latex_output_ && ! is_class_)
    {
        // Then make a file documentation block
        cout << "/** @file \"" << filename_ << "\"\n  ";

      cout_ingroup();

      cout << "*/\n";
    }
  }
  cout << "/*";
  if(latex_output_ && !is_class_)
  {
    cout_ingroup();
    cout << "\n  ";
  }
  if(is_setter_ || is_getter_)
  {
    cout << "* @var " << cfuncname_ << "\n  ";
    string temp = (is_setter_ ? "Setter" : "Getter");
    cout << "* @par " << temp << " is implemented\n  *";
  }
  else
  {
    // specify the @fn part
    cout << "* @fn ";
    print_pure_function_synopsis();

    // specify the @brief part
    cout << "\n  * @brief ";
  }
  cout_docuheader(cfuncname_);
  cout << "*\n  ";

  // specify the @details part

  // standard body definitions
  cout_docubody();

  // parameters
  if(!param_list_.empty() && !is_getter_ && !is_setter_)
  {
    cout << "*\n  ";
    write_docu_list(param_list_, "@param", cscan_.param_list_);
  }

  // return values
  if(!return_list_.empty() && !is_constructor && !is_getter_ && !is_setter_)
  {
    cout << "*\n  ";
    write_docu_list(return_list_, "@retval", cscan_.return_list_);
  }

  // required fields
  write_docu_listmap(required_list_, "@par Required fields of ", cscan_.field_docu_);

  // optional fields
  write_docu_listmap(optional_list_, "@par Optional fields of ", cscan_.field_docu_);

  // return fields
  write_docu_listmap(retval_list_, "@par Generated fields of ", cscan_.field_docu_);
  clear_lists();

  // extra docu fields
  cout_docuextra();
  if( new_syntax_ )
  {
    cout << "* @synupdate Syntax needs to be updated! \n  ";
  }
  cout << "*/\n";
  if(!is_method)
    is_first_function_ = false;

  is_setter_ = false; is_getter_ = false;
  cfuncname_.clear();
}

void usage()
{
  cout
    << "Usage: mtoc mfile [latex_output] [conf]\n"
    << "       mtoc --help\n\n"
    << "mtoc parses a Matlab m-file 'mfile'.\n"
    << "If 'latex_output' is set to 1, the output is optimized for latex output.\n"
    << "The default is 0.\n"
    << "A configuration file needs to be given or it must exist in\n"
    << "'./doxygen/mtoc.conf'." << endl;
}

// main routine
int main(int argc, char ** argv)
{
  istream * fcin;
  std::ifstream fin;
  string filename;
  if(argc >= 2)
  {
    if (string(argv[1]) == string("--help"))
    {
      usage();
      exit(0);
    }
    try
    {
      fin.open(argv[1]);
      fcin = &fin;
      filename = argv[1];
    }
    catch (std::ifstream::failure e)
    {
      cout << "Exception opening/reading file";
      usage();
      exit(-1);
    }
  }
  else
  {
    fcin = &cin;
    filename = "stdin";
  }

  bool latex_output = false;
  if(argc == 3)
  {
    latex_output = (strncmp(argv[2],"1",1)==0) ? true : false;
  }

  std::string conffilename;
  if(argc == 4)
  {
    conffilename = std::string(argv[3]);
  }

  char buf[1000];
  char * dummy = getcwd(buf, 1000);
  dummy = 0;

  string::size_type found = 0;
  string cwd(buf);
  found = filename.find(cwd);
  if(found!=string::npos)
  {
    filename = filename.substr(cwd.size()+1);
  }
  MFileScanner scanner(*fcin, filename, conffilename, latex_output);
  scanner.execute();
  return 0;
}

/* vim: set et sw=2 ft=ragel foldmethod=marker: */

