#ifndef _FEM_TOKEN_HH_
#define _FEM_TOKEN_HH_
//#pragma interface

#include <cctype>
#include <cstdlib>
#include "String.hh"
#include "Stream.hh"
#include <vector>
#include <set>

class Token
{
  
  Stream  sin;
  int     pos_prev;
  
public:
  
  enum token_value {
    none=01,function=02,name=04,option=010,floating=020,integer=040,
    phrase=0100,
    number=floating|integer,
    word=floating|integer|name|phrase
  } tvalue, tbuf, curr_tok;

  bool is_split;
  
  int ibuf;

  double fvalue; int dvalue;

  String sbuf;                      // get
  String svalue, undef;   // token

  std::vector<String> split_vector;


  Token();
  //  Token(char* str);
  Token(String str);

  void split();
  
  void attach(char* str)
    {
      sin.attach(str);
    }
  
  void attach(String& str)
    {
      sin.attach(str);
    }

  //  void csv(String str = s_undef, char field_separator = c_comma);   // transform csv remainder to token friendly
  void csv_tok(char field_separator);   // transform csv remainder to token friendly
  void xml();                                                       // transform xml remainder to std token
  
  token_value get_token();
  token_value next();
  token_value next(String& st){
    next();
    st = svalue;
    return tvalue;
  }
  int         ntoken();           // number of tokens left
  String      get_substring(int,int,             bool preserve_quote = false);
  String      get_substring(std::set<int>& v_int,bool preserve_quote = false);
  String      get_substring(std::set<int>& v_int,bool preserve_quote, bool& is_answ);
  String      del_substring(std::set<int>& v_int,bool preserve_quote = false);
  
  token_value value(){ return tvalue; }
  token_value last(String& str){
    str = svalue;
    return tvalue;
  }
  
  void reset(){                    // starts from begining
    sin.seekg(0);
  }
  
  void putback(){                  // put back token, useful for while
     sin.seekg(pos_prev);
  }
  
  String back(){                   // back substr after next command
    return sin.back();
  }

  String front(){                   // front
    return sin.front(pos_prev);
  }

  String left(){                   // back substr before next command
    return sin.back(pos_prev);
  }

  String buffer(){
    return sin.buffer();
  }
  
  String pvalue(){
    char ch= c_verbatim;
    int beg = 0,end = svalue.size();
    if(svalue[0] == ch)
      beg = 1;
    if(svalue[end-1] == ch)
      end = end-1;
    return svalue.substr(beg,end-beg);
  }
  
  token_value get(int i1,int i2); // fill sbuf with words from i1 to i2
  
  int isvalid(char ch)
    {
      if(isalnum(ch) || ch == '-' || ch == '*' ||
         ch == '~'   || ch == '@' || ch == '.' ||
         ch == ':'   || ch == '&' || ch == '+' ||
         ch == '('   || ch == ')' || ch == '/' ||
         ch == '#'   || ch == '$' || ch == '"' ||
         ch == '_'   || ch == '%' || ch == '!' ||
         ch == '<'   || ch == '>' || ch == '?' ||
         ch == '|'   || ch == '[' || ch == ']' ||
	 ch == ';'   || ch == '{' || ch == '}' ||
	 ch == '('   || ch == ')' || ch == '^' )
        return 1;
      else
        return 0;
    }
  
  int isend(char ch)
    {
      if(ch== EOF || ch == '\n'  ||  ch == '\0' )
        return 1;
      else
        return 0;
    }

  int iseat(char ch)
    {
      if(ch == ' ' || ch == '\t' || ch == '\r' || ch == '\v' ||
         ch == '\f'|| ch == ',')
        return 1;
      else
        return 0;
    }

  int issimple(char ch)
    {
      if(ch == '=')
        return 1;
      else
        return 0;
    }  

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  template< class T>  bool getkey  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  template< class T>  bool getkey(String key, T& value){
    int pos;
    return getkey(key,value,pos);
  }
  
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  template< class T> bool getkey  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  template< class T> bool getkey(String key, T& value, int& pos){
    if( ! is_split)
      split();
    
    String optkey = "-"+key, valstr;

    for(int i=0; i<split_vector.size()-1; i++){
      pos = i;
      if(split_vector[i].contains(optkey,0)){
        valstr=split_vector[i+1];
        if( valstr.put(value)){
          std::cout<<">**+++"<<value<<std::endl;
          return true;
        }
      }
    }
    std::cout<<"**---"<<std::endl;
    return false;
  }

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  bool getkey  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  bool getkey(String key){
    if( ! is_split)
      split();
    
    String optkey = "-"+key, valstr;
    std::cout<<key<<":"<<optkey<<"|<";
  
    for(int i=0; i<split_vector.size(); i++){
      if(split_vector[i].contains(optkey,0)){
        std::cout<<">**+++"<<std::endl;
          return true;
      }
    }
    std::cout<<"**---"<<std::endl;
    return false;
  }


};



#endif
