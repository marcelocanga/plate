#ifndef _SOLVER_FEM_
#define _SOLVER_FEM_

#include<string>
#include<vector>
#include<map>
#include<set>
#include "Array.hh"

class Solver {
  std::vector<double> rhs;
  std::vector<double> lhs;
public:
  void assemble();
  void solve();
};

#endif
