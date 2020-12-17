#include <iostream>
#include <fstream>
#include "Report.hh"
#include "Essential.hh"
#include "Point.hh"
#include "Solver.hh"
#include "Plate.hh"
#include "Support.hh"

Report rep;

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//
//               -----  Report::Report  -----
//
//
//
// C: Marcelo Canga. Dec 2020
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Report::Report()
{
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Report::set_fos  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Report::set_fos(std::string fin)
{
  fin.replace(fin.find(".inp"),fin.length(),".rep");
  fos.open(fin);
  fos << std::fixed << std::showpoint<< std::setprecision(6) <<std::setfill('0');

}

void Report::fos_close()
{
  fos.close();
}

std::string pad(int pad, std::string str)
{
  std::ostringstream ss;
  ss << std::left << std::setfill(' ') << std::setw(pad) << str;
  return ss.str();
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  Report& Report::operator<<(std::ostream& omanip  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Report& Report::operator<<(std::ostream& omanip(std::ostream &)){
      omanip(std::clog);
      omanip(fos);
  return *this; 
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Report::summary  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Report::summary()
{
  char po_type[] = {'n','e','n'};
  std::string str;
  
  Solver* so_pt = Solver::current();
//
//*********************************************      
//
//
//*********************************************      points
//
  rep<<" Points "<<std::endl;
  rep<<"  Name    Type   X     Y     Z     W     Rx     Ry"<<std::endl;
  rep<<"------------------------------------------------------"<<std::endl;
  for(auto const& [first,second] : Point::point_m){
    if(second->type == Point::none_t)continue;

    str = pad(10,second->name);
    rep<<
      str<<" "     <<
      po_type[second->type]<<" "<<
      second->coor(0)<<" "<<
      second->coor(1)<<" ";
    
    if(second->type == Point::node_t){
      rep<<
	"   -      "<<
	so_pt->guess(second->dof_loc) << " "<<
	so_pt->guess(second->dof_loc+1) <<
	std::endl;
      }
    else{
      rep<<
	so_pt->guess(second->dof_loc) << " "<<
	"   -      "<<
	"   -      "<<
	std::endl;
    }
    
  }
//
//*********************************************      element shear
//
  rep<<std::endl;
  rep<<std::endl;
  rep<<" Bending Moment/Shear  "<<std::endl;
  rep<<"    Element     BX     BY     BXY      SX     SY        "<<std::endl;
  rep<<"  ------------------------------------------------------"<<std::endl;

  so_pt->lhs = d_zero;
  so_pt->rhs = d_zero;
  
  for(auto const& [first,second] : Plate::plate_m){

    second->compute_stress();
    second->assemble();
    
    str = pad(10,first);
    
    for(int ip=0; ip<3; ip++){
    if(ip == 0){
    rep<<
      str          <<"   ";
    }else
      rep<<"            ";
    rep<<
      second->bending(0,ip) <<"   "<<
      second->bending(1,ip) <<"   "<<
      second->bending(2,ip) <<"   ";
    if(ip == 0){
    rep<<
      second->shear(0)      <<"   "<<
      second->shear(1);
    }
    rep << std::endl;
    }
  }
//
//*********************************************      reaction forces, moments
//
  rep<<std::endl;
  rep<<std::endl;
  prod(so_pt->lhs,so_pt->guess,so_pt->rhs);
  rep<<" Reaction Forces/Moments  "<<std::endl;
  rep<<"   Element Side   R    MX      MY    MX   MY  "<<std::endl;
  rep<<"  --------------------------------------------"<<std::endl;

  for(auto su_pt : Support::support_v){
    AInt index;
    su_pt->el_pt->get_index(su_pt->eside-1,0,Plate::support_t,index);
    str = pad(10,su_pt->el_pt->name);

    rep<<
      str<<" "<<
      so_pt->rhs(index(2))<<"  "<<
      so_pt->rhs(index(0))<<"  "<<
      so_pt->rhs(index(1))<<"  "<<
      so_pt->rhs(index(3))<<"  "<<
      so_pt->rhs(index(4))<<"  "<<
      std::endl;
  }

  //
//*********************************************      fiber rotation
//
  rep<<std::endl;
  rep<<std::endl;
  rep<<" Rotation/Deflection Gradient  "<<std::endl;
  rep<<"    Element     RX     RY         WX     WY        "<<std::endl;
  rep<<"  -------------------------------------"<<std::endl;

  for(auto const& [first,second] : Plate::plate_m){

    second->compute_stress();
    str = pad(10,first);
    
    for(int ip=0; ip<3; ip++){
    if(ip == 0){
    rep<<
      str          <<"   ";
    }else
      rep<<"             ";
    rep<<
      second->rotation(0,ip) <<"   "<<
      second->rotation(1,ip) <<"   "<<
      second->d_deflec(0,ip) <<"   "<<
      second->d_deflec(1,ip) <<"   "<<
      std::endl;
    }
  }

  //
//*********************************************      curvature
//
  rep<<std::endl;
  rep<<std::endl;
  rep<<" Curvature  "<<std::endl;
  rep<<"    Element     KX     KY     KXY        "<<std::endl;
  rep<<"  -------------------------------------"<<std::endl;

  for(auto const& [first,second] : Plate::plate_m){

    str = pad(10,first);
    second->compute_stress();
    
    for(int ip=0; ip<3; ip++){
    if(ip == 0){
    rep<<
      str          <<"   ";
    }else
      rep<<"            ";
    rep<<
      second->curvature(0,ip) <<"   "<<
      second->curvature(1,ip) <<"   "<<
      second->curvature(2,ip) <<"   "<<
      std::endl;
    }
  }

}
