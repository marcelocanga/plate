#ifndef _SOLVER_FEM_
#define _SOLVER_FEM_

#include<string>
#include<vector>
#include<map>
#include<set>
#include "Array.hh"

class Solver {

protected:
  int ndof;
  double def_poisson, def_thickness, def_young;
  Solver* current_solver;
  std::vector<double> rhs;
  std::vector<double> lhs;
  
public:

  friend class Point;
  friend class Plate;
  friend class Load;

  Solver();
  Solver* current(){ return current_solver;}
  void parse_input(std::string);
  
  void assemble();
  void solve();
  
};

#endif
