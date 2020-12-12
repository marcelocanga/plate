#ifndef __STRING_MEC_HH__
#define __STRING_MEC_HH__
//#pragma interface

#include <cmath>
#include <iostream> 
#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <complex>
#include "Factory.hh"
#include "Essential.hh"

class String;
class Regular;


class String:public std::string{
  
public:
  
  static bool is_initialize;

  enum Pattern { rxint, rxdouble, rxwhite, rxseparator };
  
  String():std::string(){}
  
  template <class T>
  String(T ch):std::string(ch){}

  String(char ch){assign(1,ch);}

  template <class T>
  String(T ch, int len):std::string(ch,len){}

  char operator [] (int ind){
    return std::string::operator[](ind);
  }

  //  void form(char* fmt ...);
  // int  scan(char* fmt ...);
  void form(std::string fmt ...);
  int  scan(std::string fmt ...);
  void serial ();

//
//************************************************      back
//
  String  back(const char ch){
    int index = find(ch);
    if(index <0 )
      return s_zero;
    else
      return substr(index+1,size()-index);
  }

  String  rback(const char ch){
    int index = rfind(ch);
    if(index <0 )
      return s_zero;
    else
      return substr(index+1,size()-index);
  }
  
  String  back(const char* ch){
    int index = find(ch);
    if(index <0 )
      return s_zero;
    else
      return substr(index+strlen(ch),size()-index);
  }

  String rback(const char* ch,int siz = 1){
    int index = rfind(ch);
    if(index <0 )
      return s_zero;
    else
      return substr(index+siz,size());
  }
  
  String rback(String str){
    int index = rfind(str);
    if(index <0 )
      return s_zero;
    else
      return substr(index+str.size(),size());
  }

  String  back(int pos){
    if(pos > 0)
      return substr(pos,size());
    else
      return substr(0,size());
  }
//
//*********************************************      concat phrase to the string
//
  void   concat(String phrase);
  String pvalue();          // returns the value of a phrase
  void   pwrap();           // warps the string into a phrase
//
//************************************************      at
//
//  String at(int i1, int i2){
//    int len = size();
//    i1 = std::min(i1, len);
//    i2 = std::min(i1+i2, len);
//    return substr(i1, i2);
//  }
//
//************************************************      front
//
  template <class T>
  String  front(T ch){
    int index = find(ch);
    if(index < 0)
      return *this;
    else if(index == 0)
      return s_zero;
    else
      return  substr(0,index);
  }
  template <class T>
  String  rfront(T ch){
    int index = rfind(ch);
    if(index < 0)
      return *this;
    else if(index == 0)
      return s_zero;
    else
      return  substr(0,index);
  }
//
//************************************************    cast   
//
  operator  char*() const{
    return (char*)c_str();
  }
  operator const char*() const{
    return (const char*)c_str();
  }
  const char* str() const{
    return c_str();
  }
//
//************************************************      contains
//
  template< class T>
  bool  contains(T ch){
    int pos = find(ch);
    if(pos >= 0)
      return true;
    else
      return false;
  }
  
  template< class T>
  bool  contains(T ch, int start){
    int pos = find(ch);
    if(pos == start)
      return true;
    else
      return false;
  }
//
//************************************************      eatwhite, leading and trailing
//
  void eatwhite(){
    int init;
    String res;
//.............................................      all whites
//
    if( (init = find_first_not_of(s_white)) < 0 )
      assign(s_zero);
    else{
      int end  = find_last_not_of(s_white);
      assign(substr(init,end-init+1));
    }
  }
//
//************************************************      eatquote
//
  void eatchar(char eatch){
    String::iterator it;
    String newstring;
    for(it = this->begin(); it != this->end(); it++)
      if(*it != eatch)
        newstring += *it;
    assign(newstring);
  }
//
//*********************************************      case, instead of uoc::case
//
  void macro_case();
  void upper_case();
  void lower_case();
  void capital_case();
  void quote_case();
  void token_case();
  void sort();
//
//************************************************      eatquote
//
  void find_replace(char before,char after);
  void find_replace(String& before, String& after);
//
//************************************************      match
//
  bool  match(const char* ch){
    if(compare(ch) == 0)
      return true;
    else
      return false;
  }

  bool match(Regular& reg);
  bool  match(Pattern rx);
//
//************************************************      overlap
//
  static String overlap(String str1,String str2){
    int i1 = str1.length() - str2.length();
    if( i1 > 0 )
      return str1.substr(0,i1)+str2;
    else
      return str2;
  }
//
//************************************************      transform string to object T
//
  template< class T>
  bool put(T& obj){
    return obj.put(*this);
  }
    
  
  bool put(String& str){
    str= c_str();
    return true;
  }

  bool put(bool& val){

    if(*this == "1" || *this == s_true || *this == s_yes)
      val = true;
    else if(*this == "0" || *this == s_false || *this == s_no )
      val = false;
    else
      return false;

    return true;
  }

  bool put(int& val){
    char *ptr;
    int ival = strtol(c_str(),&ptr,10);
    if( *ptr == '\0'){
      val = ival;
      return true;
    }
    else
      return false;
  }

  bool put(long int& val){
    char *ptr;
    long int ival = strtol(c_str(),&ptr,10);
    if( *ptr == '\0'){
      val = ival;
      return true;
    }
    else
      return false;
  }
  
  bool put(float& val){
    char *ptr;
    double ival = strtof(c_str(),&ptr);
    if( *ptr == '\0'){
      val = ival;
      return true;
    }
    else
      return false;
  }

  bool put(double& val){
    char *ptr;
    double ival = strtod(c_str(),&ptr);
    if( *ptr == '\0'){
      val = ival;
      return true;
    }
    else
      return false;
  }

  bool put( std::complex<double> & val) {
    return false;
  }

//
//************************************************      transform  object T to string
//
  template< class T>
  bool get(T& obj){
    return obj.get(*this);
  }
    
  bool get(String& str){
    form("%s", str.c_str());
    return true;
  }
  
  bool get(bool val){
    if(val)
      form("%s",s_true);
    else
      form("%s",s_false);
    return true;
  }
  
  bool get(const char *st){
    form("%s", st);
    return true;
  }

  bool get(int val){
    form("%i",val);
    return true;
  }

  bool get(long int val){
    form("%ld",val);
    return true;
  }
  
  bool get(double val){
    form("%g",val);
    return true;
  }

  bool get(float val){
    form("%g",val);
    return true;
  }

  bool get(std::complex<double> val){
    form("%g %g i",val.real(),val.imag());
    return true;
  }

//
//************************************************      prepend
//
  void prepend(String& prep){
    assign(prep+*this);
  }
//
//************************************************      init
//
  static bool initialize(){
    Factory<String>::initialize();
    return true;
  }
//
//************************************************      done
//
  
};

#endif
