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
  double area, d_area, pressure;
  int nnode,nedge,nidof,nshear,eldim,nedof;
  double poisson, thickness, young;

  std::vector<Point*> point_v;
  std::vector<int>    edof_loc;

  static int ip;
  static double wgt,xr,xs;
  static double wg[3], xg[2][3];
  static ADouble shape,shape_h;
  static MDouble d_shape, d_shape_h;
  static MDouble b_grad, w_grad;
  static MDouble constitutive_b;
  static MDouble constitutive_s;
  static ADouble fint;
  static MDouble stiff;

public:
  friend class Point;
  friend class Solver;
  friend class Load;

  enum IndexType { support_t, moment_t, force_t};

  static std::map<std::string, Plate*> plate_m;
  
  Plate();
  void init();
  void SamplePoint(int);
 
  void add_edge();
  void count_dof(int& icount);
  
  void assemble();
  void potential();

  void compute_constitutive();
  void Grad(int);
  void Shape();

  void Fint();
  void Stiffness();

  void get_index(int side,int ind, enum IndexType, AInt& index);
};


#endif
