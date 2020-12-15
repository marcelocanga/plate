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
  static Solver* current_solver;
  MDouble lhs; 
  ADouble rhs;
public:

  friend class Point;
  friend class Plate;
  friend class Load;
  friend class Support;

  Solver();
  static Solver* current(){ return current_solver;}
  void parse_input(std::string);
  
  void run();
  void assemble();
  void solve(MDouble& mat, ADouble& vec);

  void add_rhs(AInt&, double);
  void set_ans(AInt&, double);
  
};

#endif
