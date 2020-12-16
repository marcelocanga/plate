#ifndef _LOAD_FEM_
#define _LOAD_FEM_

#include<string>
#include<vector>
#include<map>
#include<set>
#include "Array.hh"
#include "Plate.hh"
#include "Point.hh"
class Load {
  
protected:

  Plate* el_pt;
  Point*  po_pt;
  int eside, gdir;
  double value;
  
public:
  
  friend class Point;
  friend class Plate;
  friend class Solver;
  
  enum Type {moment_a,force_a,line_moment_a,line_force_a} type;
  static std::vector<Load*> load_v;
  
  Load(enum Type _type) : type(_type) {}

  void assemble();
};


#endif
