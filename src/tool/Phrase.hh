#ifndef _ME_PHRASE_HH_
#define _ME_PHRASE_HH_

#include "String.hh"

class Phrase : public String {

public:
  Phrase():String(){}
  
  template <class T>
  Phrase(T ch):String(ch){}

  Phrase(char ch){assign(1,ch);}

  template <class T>
  Phrase(T ch, int len):String(ch,len){}

  char operator [] (int ind){
    return String::operator[](ind);
  }

  bool put(String&);
  bool get(String&);
  String concat();
  
};

#endif
