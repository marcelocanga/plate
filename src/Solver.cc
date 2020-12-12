#include "Plate.hh"
#include "Load.hh"
#include "Solver.hh"

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//
//               -----  Solver::Solver  -----
//

//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Solver::Solver()
{
  current_solver = this;
  poisson        = 0.3;
  young          = 2.e6;
  thickness       = 0.1;
}



//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Solver::parse_input  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Solver::parse_input(std::string file_name){
  
  enum { point_a,plate_a,force_a,moment_a,support_a,solve_a,end_a };
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
    else if (line.find("*plat")   != std::string::npos) token = plate_a;
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
//.............................................      plate
//
    case plate_a:{
      std::vector<std::string> point(3);
      Plate* el_pt = new Plate;

      iss >> el_pt->name;
      iss >> point[0]>>point[1]>>point[2];

      Plate::plate_m[el_pt->name] = el_pt;
      
      std::cout<<"Plate:"<<el_pt->name<<","<<point[0]<<","<<point[1]<<","<<point[2]<<std::endl;
    }
      break;
//.............................................      shear force
//
    case force_a:{
      Load *lo_pt = new Load(Load::force_a);
      iss >> lo_pt->ename >>lo_pt->eside >> lo_pt->value;
      Load::load_v.push_back(lo_pt);

      std::cout<<"Shear Force. Plate:"<<lo_pt->ename <<",side:"<<lo_pt->eside <<",value:"<<lo_pt->value <<std::endl;
    }
      break;
//.............................................      moment
//
    case moment_a:{
      Load *lo_pt = new Load(Load::moment_a);
      iss >> lo_pt->ename >>lo_pt->eside >> lo_pt->value;
      Load::load_v.push_back(lo_pt);
      std::cout<<"Bending Moment. Plate:"<<lo_pt->ename <<"side:"<<lo_pt->eside <<",value:"<<lo_pt->value<<std::endl;
    }
      break;
//.............................................      support
//
    case support_a:{
      Load *lo_pt = new Load(Load::support_a);
      iss >> lo_pt->ename >>lo_pt->eside;
      Load::load_v.push_back(lo_pt);
      std::cout<<"Fix Support. Plate:"<<lo_pt->ename <<"side:"<<lo_pt->eside <<std::endl;
    }
      break;
//.............................................      solver
//
    case solve_a:{
      for(auto const& [first,second] : Plate::plate_m) second->add_point();
      assemble();
      solve();
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


void Solver::assemble()
{
  //  for(auto const& [first,second] : Plate::element_m) second->build_index();
  for(auto const& [first,second] : Plate::plate_m) second->assemble();
  for(auto const&  first         : Load::load_v)       first->assemble();
}     

void Solver::solve()
{}
