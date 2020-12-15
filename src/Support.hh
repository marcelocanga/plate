#ifndef _SUPPORT_FEM_
#define _SUPPORT_FEM_

#include<string>
#include<vector>
#include<map>
#include<set>
#include "Array.hh"
#include "Plate.hh"

class Support {
  
protected:

  Plate* el_pt;
  int eside;
  
public:
  
  friend class Point;
  friend class Plate;
  friend class Solver;
  
  static std::vector<Support*> support_v;
  
  void assemble();

};


#endif
