#ifndef _LOAD_FEM_
#define _LOAD_FEM_

#include<string>
#include<vector>
#include<map>
#include<set>
#include "Array.hh"

class Load {
public:
  enum Type {moment_a,force_a,support_a} type;
  std::string ename;
  int eside;
  double value;
  static std::vector<Load*> load_v;
  Load(enum Type _type) : type(_type) {}
  void assemble();
  void potential();
};


#endif
