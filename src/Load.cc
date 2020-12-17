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
// C: Marcelo Canga. Dec 2020
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


void Load::assemble()
{
  bool is_index = false;
  double fac,length;
  AInt index;
  
  int nd1     =  eside-1;
  int nd2     =  eside %3;

  switch(type){
  case moment_a:
    is_index = po_pt->get_index(gdir-1,index);
    fac      = value;
    break;

  case line_moment_a:
    length   = norm_2(el_pt->point_v[nd1]->coor,el_pt->point_v[nd2]->coor);
    is_index = el_pt->get_index(eside-1,gdir-1,Plate::moment_t,index);
    fac = 0.5 * value * length;
    break;

  case force_a:
  case line_force_a:
    length   = norm_2(el_pt->point_v[nd1]->coor,el_pt->point_v[nd2]->coor);
    is_index = el_pt->get_index(eside-1,0,Plate::force_t,index);
    fac = value ;
    if(type == line_force_a) fac = value * length;
    break;
  }
  
  if(is_index) Solver::current()->add_rhs(index,fac);

}
