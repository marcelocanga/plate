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
  
  static int ninteg;
  double area, d_area, pressure;
  int nnode,nedge,nidof,nshear,eldim,nedof;
  double poisson, thickness, young;

  std::vector<Point*> point_v;
  std::vector<int>    edof_loc;

  static double wgt,xr,xs;
  //  static ADouble wg;
  //  static MDouble xg;
  static std::vector<double> xg_v[2];
  static std::vector<double> wg_v;
  static ADouble shape,shape_h;
  static MDouble d_shape, d_shape_h;
  static MDouble b_grad, w_grad;
  static MDouble constitutive_b;
  static double  constitutive_s;
  static ADouble fint;
  static MDouble stiff;

  static MDouble bending,rotation,d_deflec,curvature;
  static ADouble shear;

public:
  friend class Point;
  friend class Solver;
  friend class Load;
  friend class Report;

  enum   IndexType   { support_t, moment_t, force_t};
  enum   ShapeOrigin { corner_t = 1, edge_t = 2 };
  static ShapeOrigin shape_origin;
  
  static std::map<std::string, Plate*> plate_m;
  
  Plate();
  void init();
  static void setup();
  static void integration_point();
  
  void SamplePoint(int);
 
  void add_edge();
  void count_dof(int& icount);
  
  void assemble();
  void potential();

  void compute_constitutive();
  void Grad(int);
  void Shape();
  void ShapeCorner();
  void ShapeEdge();

  void Fint();
  void Stiffness();

  bool get_index(int side,int ind, enum IndexType, AInt& index);

  void compute_stress();
  void Stress(int);
};


#endif
