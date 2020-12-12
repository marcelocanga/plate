#include<iostream>
#include<fstream>
#include<sstream>
#include "Plate.hh"

extern "C"{
  int test_();
  int test2_();
}

extern "C"{
#include "Blas.h"
#include "Lapack.h"
}


std::map<std::string, Point*>    Point::point_m;
std::map<std::string, Element*>  Element::element_m;
std::set<Point*,LtPoint>         Point::u_point_s;
std::vector<Load*>               Load::load_v;

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void parse_input  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void parse_input(std::string file_name){
  
  enum { point_a,element_a,force_a,moment_a,support_a,solve_a,end_a };
  int token = end_a;

  std::string line;
  std::ifstream input;
  
  input.open(file_name);
//
//*********************************************      switch
//
  while(std::getline(input,line)){
    std::istringstream iss(line);
    bool is_continue = true;

    if      (line.find("$")       != std::string::npos) continue;
    else if (line.find("*poin")   != std::string::npos) token = point_a;
    else if (line.find("*elem")   != std::string::npos) token = element_a;
    else if (line.find("*forc")   != std::string::npos) token = force_a;
    else if (line.find("*mome")   != std::string::npos) token = moment_a;
    else if (line.find("*supp")   != std::string::npos) token = support_a;
    else if (line.find("*solv")   != std::string::npos) token = solve_a;    
    else if (line.find("*end")    != std::string::npos) break;
    else is_continue = false;

    if(is_continue) continue; 
    
    switch(token){
//.............................................      point
//
    case point_a:{
      Point* po_pt = new Point;
      iss >> po_pt->name;
      iss >> po_pt->coor(0) >> po_pt->coor(1) >> po_pt->coor(2) ;
      Point::point_m[po_pt->name] = po_pt;
      std::cout<<"Point:"<<po_pt->name<<":"<<po_pt->coor(0)<<","<<po_pt->coor(1)<<","<<po_pt->coor(2)<<std::endl;
    }
      break;
//.............................................      element
//
    case element_a:{
      std::vector<std::string> point(3);
      Element* el_pt = new Element;

      iss >> el_pt->name;
      iss >> point[0]>>point[1]>>point[2];

      Element::element_m[el_pt->name] = el_pt;
      
      std::cout<<"Element:"<<el_pt->name<<","<<point[0]<<","<<point[1]<<","<<point[2]<<std::endl;
    }
      break;
//.............................................      shear force
//
    case force_a:{
      Load *lo_pt = new Load(Load::force_a);
      iss >> lo_pt->ename >>lo_pt->eside >> lo_pt->value;
      Load::load_v.push_back(lo_pt);

      std::cout<<"Shear Force. Element:"<<lo_pt->ename <<",side:"<<lo_pt->eside <<",value:"<<lo_pt->value <<std::endl;
    }
      break;
//.............................................      moment
//
    case moment_a:{
      Load *lo_pt = new Load(Load::moment_a);
      iss >> lo_pt->ename >>lo_pt->eside >> lo_pt->value;
      Load::load_v.push_back(lo_pt);
      std::cout<<"Bending Moment. Element:"<<lo_pt->ename <<"side:"<<lo_pt->eside <<",value:"<<lo_pt->value<<std::endl;
    }
      break;
//.............................................      support
//
    case support_a:{
      Load *lo_pt = new Load(Load::support_a);
      iss >> lo_pt->ename >>lo_pt->eside;
      Load::load_v.push_back(lo_pt);
      std::cout<<"Fix Support. Element:"<<lo_pt->ename <<"side:"<<lo_pt->eside <<std::endl;
    }
      break;
//.............................................      solver
//
    case solve_a:{
      Solution sol;
      for(auto const& [first,second] : Element::element_m) second->add_point();
      sol.assemble();
      sol.solve();
    }
      break;
//
//*********************************************      end loop,switch
//
    }
  }
//
//*********************************************      done
//
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  int main  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

int main(int argc, char *argv[])
{
  std::string input_file;
  
  if(argc > 1) input_file = argv[1];
  else         input_file = "plate.inp";

  parse_input(input_file);
  
}
