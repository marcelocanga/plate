#ifndef _PLATE_FEM_
#define _PLATE_FEM_

#include<string>
#include<vector>
#include<map>
#include<set>
#include "Array.hh"

class Plate{
public:
  int ninteg;
  double young, poisson, tickness, area;
  
  std::string name;
  std::vector<Point*> point_v;
  std::vector<int>    index_v;
  static std::map<std::string, Element*> element_m;

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
