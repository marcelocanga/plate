
#ifndef __FACTORY_MEC_HH__
#define __FACTORY_MEC_HH__
//#pragma interface

#include <iostream>
#include <fstream>
#include <map>
#include "Diag.hh"

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//
//               -----  class TypeFactory  -----
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
class _TypeFactory{
protected:
  static int next_id;
};

template <class Derived>
class TypeFactory : private _TypeFactory{
  static int id;
public:
  static int set_id(){
    id = id ? id : (id = ++next_id);
    return id;
  }
  static int get_id(){
    return id;
  }
};

template <class Derived> int TypeFactory<Derived>::id = 0;


//
//************************************************   one factory per basic type
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//
//               -----  class Factory  -----
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

template <class Basic>
class Factory{
  static int type_id;
  static int derived_size;
  static bool is_polymorphic_;
public:
//
//************************************************      map type to new()
//
  static std::map<int, Basic* (*) () >&
  registry(){
    static std::map<int, Basic* (*) ()>* pmap =
      new  std::map<int, Basic* (*) ()>;
    return *pmap;
  }
  static bool is_polymorphic(){
    return is_polymorphic_;
  }
//
//************************************************      function for map
//
  template <class Derived>
  static Basic* create(){
    derived_size = sizeof(Derived);
    return new Derived();
  }
//
//************************************************      initialization
//
  static bool initialize() {
    is_polymorphic_ = false;
    return true;
  }
  
  template <class Derived>
  static bool initialize() {
    is_polymorphic_ = true;
    if(is_polymorphic_){
      int index = TypeFactory<Derived>::set_id();
      Factory::registry()[index] = Factory::template create<Derived>;
    }
    return true;
  }
//
//************************************************      done
//
};

template <class Basic> int  Factory<Basic>::type_id;
template <class Basic> int  Factory<Basic>::derived_size = sizeof(Basic);
template <class Basic> bool Factory<Basic>::is_polymorphic_(false);

#endif
