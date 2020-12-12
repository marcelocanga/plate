#ifndef _PLATE_FEM_
#define _PLATE_FEM_

#include<string>
#include<vector>
#include<map>
#include<set>

#include "Array.hh"
#include "Point.hh"

class Plate{

protected:
  
  int ninteg;
  double area;
  
  std::string name;
  std::vector<Point*> point_v;
  std::vector<int>    index_v;

public:
  friend class Point;
  friend class Solver;
  friend class Load;
  
  static std::map<std::string, Plate*> plate_m;
  
  void init();
  void SamplePoint(int);
 
  void add_point();
  void assemble();
  void potential();

  void compute_constitutive();
  void Grad(int);
  void Shape();
  void ShapeH();

  void Force();
  void Fint();
  void Stretch();
  void Stiffness();
};


#endif
