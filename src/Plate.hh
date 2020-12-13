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
  
  std::string name;
  
  int ninteg;
  double area;
  int nnode,nedge,nidof,eldim,ndof;
  double poisson, thickness, young;

  std::vector<Point*> point_v;
  AInt edof_loc;

  static int ip;
  static double wgt,xr,xs;
  static double wg[3], xg[3][2];
  static ADouble shape,shape_h;
  static MDouble b_shape, b_shape_h;
  
  static ADouble fint;
  static MDouble stiff;

public:
  friend class Point;
  friend class Solver;
  friend class Load;
  
  static std::map<std::string, Plate*> plate_m;
  
  void init();
  void SamplePoint(int);
 
  void add_edge();
  void count_dof(int& icount);
  
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
