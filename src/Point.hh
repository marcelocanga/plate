#ifndef _POINT_FEM_
#define _POINT_FEM_

#include<string>
#include<vector>
#include<map>
#include<set>
#include "Array.hh"

struct LtPoint;

class Point {
  
protected:
  
  std::string name;
  Coord coor;
  int dof_loc;
  
public:
  
  friend class LtPoint;
  friend class Solver;
  friend class Plate;
  friend class Load;

  Point();

  static std::map<std::string, Point*> point_m;
  static std::set<Point*,LtPoint>      u_point_s;

  static std::string new_name();

  
};


struct LtPoint
{
  bool operator()(Point* nodp1, Point* nodp2) const {
    bool  answ = false;
    double eps = 1.e-8;

    if((int)(nodp1->coor(0)/eps)+1 <  (int)(nodp2->coor(0)/eps)){
      answ = true;
      goto done;
    }
    if((int)(nodp1->coor(0)/eps)-1 >  (int)(nodp2->coor(0)/eps)){
      answ = false;
      goto done;
    }
    if((int)(nodp1->coor(1)/eps)+1 <  (int)(nodp2->coor(1)/eps)){
      answ = true;
      goto done;
    }
    if((int)(nodp1->coor(1)/eps)-1 >  (int)(nodp2->coor(1)/eps)){
      answ = false;
      goto done;
    }
    if((int)(nodp1->coor(2)/eps)+1 <  (int)(nodp2->coor(2)/eps)){
      answ = true;
      goto done;
    }
  done:
    return answ;
  }
};



#endif
