#ifndef __MECOM_STREAM_HH__
#define __MECOM_STREAM_HH__

#include "String.hh"

class Stream{

  String str;
  int ipos, length;

public:
  
  Stream():ipos(0),length(0){}
  
  Stream(String st):str(st),ipos(0),length(str.length()){}

  Stream(char* st):str(st),ipos(0),length(str.length()){}

  Stream(char* st, int len):str(st,len),ipos(0),length(len){}

  void attach(String& st){
    ipos=0;
    str=st;
    length=str.length();
  }

  void attach(char* st){
    ipos=0;
    str=st;
    length=str.length();
  }

  void attach(char* st, int len){
    String x(st,len);
    str=x;
    length=str.length();
    ipos=0;
  }
  
  void attach(int beg, String& saf){      // attach saf after beg

    if(beg >= 0 && beg < str.length())
      str    = str.substr(0,beg) + saf;
    else
      str    += saf;
       
    length = str.length();
    ipos   = beg;
       
  }

  int get(){
    if(ipos<length){
      char ch  = str[ipos++];
      return ch;
    }
    else{
      ipos++;
      return EOF;
    }
  }

  int pos(){
    return ipos;
  }

  String buffer(){
    return str;
  }
  
  String front(){
    return front(ipos);
  }
  
  String front(int curr_pos ){
    String fro;

    if(curr_pos < str.length()){
      if(curr_pos > 0)
        fro=str.substr(0,curr_pos);
      else
        fro=str.substr(0,str.size());
    }
    else
      fro="";
    return fro; 
  }

  String back(){
    return back(ipos);
  }

  String substr(int ibeg, int isize){
    return str.substr(ibeg,isize);
  }
  
  String back(int curr_pos ){
    String aft;
    if(curr_pos < str.length()){
      if(curr_pos > 0)
        aft=str.substr(curr_pos,str.size());
      else
        aft=str.substr(0,str.size());
    }
    else
      aft="";
    return aft; 
  }
  
  int unget(){
    ipos--;
    if(ipos <0)
      return EOF;
    else
      return 1;
  }
    
  int seekg(int ind){
    if(ind < length){
      ipos = ind;
      return 1;
    }
    else
      return EOF;
  }


};


#endif
