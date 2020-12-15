#include "Plate.hh"
#include "Load.hh"
#include "Solver.hh"

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//
//               -----  void Load::assemble  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


void Load::assemble()
{
  double fac;
  AInt index;
  
  int nd1     =  eside-1;
  int nd2     =  eside %3;

  double length = norm_2(el_pt->point_v[nd1]->coor,el_pt->point_v[nd2]->coor);

  switch(type){
  case moment_a:
    el_pt->get_index(eside-1,gdir-1,Plate::moment_t,index);
    fac = 0.5 * value * length;
    break;
  case force_a:
    el_pt->get_index(eside-1,0,Plate::force_t,index);
    fac = value * length;
    break;
  }
  
  Solver::current()->add_rhs(index,fac);

}
