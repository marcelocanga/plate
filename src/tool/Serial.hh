#ifndef __SERIAL_MEC_HH__
#define __SERIAL_MEC_HH__


#include <iostream>
#include <fstream>
#include <map>
#include "Diag.hh"
#include "Factory.hh"

class Serial;
extern Serial serializer;


//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//
//   -----  void quiteserial(char action, Obj* obj, iostream& fs)  -----
//
//					    Copyright
//                                          Marcelo Canga
//					    Year 2000-2004
//                                          All Rights Reserved
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class Serial{

  int  len;
  char action;
  std::iostream* fs;

public:
  
  enum Format {binary,ascii} format;

  char get_action(){
    return action;
  }
  void length(int size){
    len += size;
  }

  std::iostream& stream(){
    return *fs;
  }
    
  
  template< class T>
   Serial& operator , (T& p){
    p.serial();
    return *this;
  }
  
  template< class T>
   Serial& operator , (T*& p){
//
//************************************************      set
//
    bool is_polymorphic = Factory<T>::is_polymorphic();
    int type_id;
//
//************************************************      switch
//
    switch(action){
//
//************************************************      read
//
    case 'r':{

   // if(p != 0) delete p;
      
      if(is_polymorphic){
        serializer<type_id;
        p =  Factory<T>::registry()[type_id]();
      }
      else
        p = new T;

      serializer<*p;
      break;
    }
//
//************************************************      write
//
    case 'w':
      if(is_polymorphic){
        type_id = p->get_serial_id();
        serializer<type_id;
      }
      serializer<*p;
      break;
//
//************************************************      size
//
    case 's':
      break;
    }
    return *this;
//
//************************************************      done
//
  }
  
  Serial& operator () (char action_, std::iostream& fs_, 
                       enum Format format_=binary){
    len = 0;
    fs = &fs_;
    action = action_;
    format = format_;
    return *this;
  }
  
  Serial& operator () (char action_, std::iostream& fs_, int& len_,
                       enum Format format_=binary){
    len = len_;
    fs = &fs_;
    action = action_;
    format = format_;
    return *this;
  }

  template< class T>
  Serial& operator < (T& p){
    *this , p;
    return *this;
  }
  
  template <class Obj>
  void raw_single(Obj& obj);
  template <class Obj>
  void raw_vector(Obj*& obj,int num);
  template <class Obj>
  void pointer(Obj *& obj, int num);
  
};


template <class Obj> inline
void Serial::raw_single(Obj& obj){
  int size;
  char ch;
  switch(format){
  case binary:{
    size = sizeof(Obj);
    switch(action){
    case 'r':{
      fs->read((char*)&obj,size);
      break;
    }
    case 'w':
      fs->write((char*)&obj,size);
      break;
    default:
      diag_mesg(diag.error,"Serialize call without valid action\n");
    }
  }
    break;
  case ascii:
    switch(action){
    case 'r':
      *fs>>obj>>ch;
      break;
    case 'w':
      *fs<<obj<<c_space;
      break;
    default:
      diag_mesg(diag.error,"Serialize call without valid action\n");
    }
    break;
  }    
}

template <class Obj> inline
void Serial::raw_vector(Obj*& obj, int num){
  int size;
  char ch;
  switch(format){
  case binary:{
    size = sizeof(Obj);
    switch(action){
    case 'r':{
      //if(obj != 0) delete[] obj;
      obj = new Obj[num];
      fs->read((char*)obj,size*num);
      break;
    }
    case 'w':
      fs->write((char*)obj,size*num);
      break;
    case 's':
      length(sizeof(Obj)*num);
      break;
    default:
      diag_mesg(diag.error,"Serialize call without valid action\n");
    }
  }
    break;
  case ascii:
    switch(action){
    case 'r':
      //if(obj) delete[] obj;
      obj = new Obj[num];
      for(int i=0; i<num; i++){
        *fs>>*(obj+i)>>ch;
      }
      break;
    case 'w':
      for(int i=0; i<num; i++)
        *fs<<*(obj+i)<<c_space;
      break;
    case 's':
      length(sizeof(Obj)*num);
      break;
    default:
      diag_mesg(diag.error,"Serialize call without valid action\n");
    }
    break;
  }    
}

template <class Obj> inline
void Serial::pointer(Obj*& obj, int num){
  
  if(action ==  'r'){
   //if(obj) delete[] obj;
    obj = new Obj[num];
  }
  
  for(int i = 0; i< num; i++)
    serializer<*(obj+i);
  
}


template <> inline
void Serial::pointer(bool*& p, int num){
  raw_vector(p,num);
}


template <> inline 
void Serial::pointer(char*& p, int num){
  raw_vector(p,num);
}

template <> inline 
void Serial::pointer(int*& p, int num){
  raw_vector(p,num);
}

template <> inline 
void Serial::pointer(long*& p, int num){
  raw_vector(p,num);
}

template <> inline 
void Serial::pointer(float*& p, int num){
  raw_vector(p,num);
}

template <> inline 
void Serial::pointer(double*& p, int num){
  raw_vector(p,num);
}



template<> inline
 Serial& Serial::operator , (bool& p){
  raw_single(p);
  return *this;
}
template<>  inline
 Serial& Serial::operator , (char& p){
  raw_single(p);
  return *this;
}
template<> inline
 Serial& Serial::operator , (int& p){
  raw_single(p);
  return *this;
}
template<> inline
 Serial& Serial::operator , (long& p){
  raw_single(p);
  return *this;
}
template<>  inline
 Serial& Serial::operator , (float& p){
  raw_single(p);
  return *this;
}
template<>  inline
 Serial& Serial::operator , (double& p){
  raw_single(p);
  return *this;
}


#endif
