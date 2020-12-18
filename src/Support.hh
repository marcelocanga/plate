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
  int eside, gdir;
  
public:
  
  friend class Point;
  friend class Plate;
  friend class Solver;
  friend class Report;

  enum Type { fix_t, line_t} type;

  Support();
  
  static std::vector<Support*> support_v;
  
  void assemble();

};


#endif
