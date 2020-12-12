#ifndef _PLATE_FEM_
#define _PLATE_FEM_

#include<string>
#include<vector>
#include<map>
#include<set>
#include "Array.hh"

struct LtPoint;

class Point {
public:
  std::string name;
  Coord coor;
  static std::map<std::string, Point*> point_m;
  static std::set<Point*,LtPoint> u_point_s;
  std::string new_name();
  
};


struct LtPoint
{
  bool operator()(Point* nodp1, Point* nodp2) const {
    bool  answ = false;
    double eps = 1.e-8;

    if((int)(nodp1->coor[0]/eps)+1 <  (int)(nodp2->coor[0]/eps)){
      answ = true;
      goto done;
    }
    if((int)(nodp1->coor[0]/eps)-1 >  (int)(nodp2->coor[0]/eps)){
      answ = false;
      goto done;
    }
    if((int)(nodp1->coor[1]/eps)+1 <  (int)(nodp2->coor[1]/eps)){
      answ = true;
      goto done;
    }
    if((int)(nodp1->coor[1]/eps)-1 >  (int)(nodp2->coor[1]/eps)){
      answ = false;
      goto done;
    }
    if((int)(nodp1->coor[2]/eps)+1 <  (int)(nodp2->coor[2]/eps)){
      answ = true;
      goto done;
    }
  done:
    return answ;
  }
};


class Element{
public:
  int ninteg;
  double young, poisson, tickness, area;
  
  std::string name;
  static std::map<std::string, Element*> element_m;
  std::vector<Point*> point_v;
  std::vector<int>    index_v;

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

class Load {
public:
  enum Type {moment_a,force_a,support_a} type;
  std::string ename;
  int eside;
  double value;
  static std::vector<Load*> load_v;
  Load(enum Type _type) : type(_type) {}
  void assemble();
  void potential();
};

class Solution {
  std::vector<double> rhs;
  std::vector<double> lhs;
public:
  void assemble();
  void solve();
};

#endif
