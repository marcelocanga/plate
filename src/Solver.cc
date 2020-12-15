#include "Plate.hh"
#include "Load.hh"
#include "Solver.hh"
#include "Support.hh"

extern "C"{      // Lapack solvers
  //#include "Blas.h"
  //#include "Lapack.h"
  int dposv_ (char*,int*,int*,double*,              int*,     double*,               int*, int*);
  int dgesv_ (int*,int*,      double*,              int*,int*,double*,               int*, int*);
}

Solver* Solver::current_solver;

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
  ndof = 0;
  current_solver = this;
  def_poisson        = 0.3;
  def_young          = 2.e6;
  def_thickness      = 0.1;
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
  
  enum { point_a,plate_a,force_a,moment_a,pressure_a,support_a,solve_a,young_a,
         poisson_a, thickness_a, end_a };
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
    else if (line.find("*pres")   != std::string::npos) token = pressure_a;
    else if (line.find("*supp")   != std::string::npos) token = support_a;
    else if (line.find("*solv")   != std::string::npos) token = solve_a;    
    else if (line.find("*youn")   != std::string::npos) token = young_a;    
    else if (line.find("*pois")   != std::string::npos) token = poisson_a;    
    else if (line.find("*thic")   != std::string::npos) token = thickness_a;    
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

      for (auto nam : point)
      {
        if (Point::point_m.count(nam))
          el_pt->point_v.push_back(Point::point_m[nam]);
        else
        {
          std::cout << "Point in plate not found:" << nam << ", plate" << el_pt->name << std::endl;
          break;
        }
      }
      Plate::plate_m[el_pt->name] = el_pt;

      el_pt->poisson = def_poisson;
      el_pt->young = def_young;
      el_pt->thickness = def_thickness;

      std::cout<<"Plate:"<<el_pt->name<<","<<point[0]<<","<<point[1]<<","<<point[2]<<std::endl;
    }
      break;
//.............................................      shear force
//
    case force_a:{
      std::string ename;
      Load *lo_pt = new Load(Load::force_a);
      iss >> ename >>lo_pt->eside >> lo_pt->value;

      if(Plate::plate_m.count(ename))
        lo_pt->el_pt = Plate::plate_m[ename];
      else{
        std::cout << "Plate in load force not found:" << ename << std::endl;
        break;
      }

      Load::load_v.push_back(lo_pt);

      std::cout<<"Shear Force. Plate:"<<ename <<", side:"<<lo_pt->eside <<", value:"<<lo_pt->value <<std::endl;
    }
      break;
//.............................................      moment
//
    case moment_a:{
      int gdir;
      std::string ename;
      Load *lo_pt = new Load(Load::moment_a);
      
      iss >> ename >>lo_pt->eside >> gdir >> lo_pt->value;
      lo_pt->gdir = gdir-1;  // index from 0

      if(Plate::plate_m.count(ename))
        lo_pt->el_pt = Plate::plate_m[ename];
      else{
        std::cout << "Plate in load moment not found:" << ename << std::endl;
        break;
      }

      Load::load_v.push_back(lo_pt);

      std::cout<<"Bending Moment Load. Plate:"<<ename <<", side:"<<lo_pt->eside <<", dir:"<<gdir<<", value:"<<lo_pt->value<<std::endl;
}
      break;
//.............................................      pressure
//
    case pressure_a:{
      std::string ename;
      double val;
      iss >> ename >> val;

      if(ename == "@"){
	for(auto const& [first,second] : Plate::plate_m){
	  second->pressure += val;
	}
	break;
      }
      
      if(Plate::plate_m.count(ename)){
	Plate::plate_m[ename]->pressure += val;
      }
      else{
        std::cout << "Plate in pressure not found:" << ename << std::endl;
        break;
      }

    }
      break;
//.............................................      support
//
    case support_a:{
      std::string ename;
      Support *su_pt = new Support;
      iss >> ename >>su_pt->eside;

      if(Plate::plate_m.count(ename))
        su_pt->el_pt = Plate::plate_m[ename];
      else{
        std::cout << "Plate in support not found:" << ename << std::endl;
        break;
      }

      Support::support_v.push_back(su_pt);

      std::cout<<"Fix Support Condition. Plate:"<<ename <<", side:"<<su_pt->eside <<std::endl;
    }
      break;
//.............................................      material properties/thickness
//
    case thickness_a:
      iss >> def_thickness;
      std::cout<<"Plate thickness:"<<def_thickness<<std::endl;
      break;
      
    case young_a:
      iss >> def_young;
      std::cout<<"Young's modulus:"<<def_young<<std::endl;
      break;

    case poisson_a:
      iss >> def_poisson;
      std::cout<<"Poisson's ratio:"<<def_poisson<<std::endl;
      break;
//.............................................      solver
//
    case solve_a:{
      run();
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
//               -----  void Solver::run  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Solver::run()
{
  ndof = 0;
  std::cout<<"Solver starts."<<std::endl;
//
//*********************************************     add nodes at element edges 
//
  for(auto const& [first,second] : Plate::plate_m) second->add_edge();
//
//*********************************************      count and map equation index
//
  for(auto const& [first,second] : Plate::plate_m) second->count_dof(ndof);
  std::cout<<"Ndof:"<<ndof<<std::endl;
//
//*********************************************      assemble stifness matrix and force vector
//
  assemble();
//
//*********************************************      solve system of equations
//
  solve(lhs,rhs);
//
//*********************************************      done
//
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Solver::assemble  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Solver::assemble()
{
  lhs.dim(ndof,ndof);
  rhs.dim(ndof);
//
//*********************************************      elements
//
  for(auto const& [first,second] : Plate::plate_m)     second->assemble();
//
//*********************************************      loads
//
  for(auto const&  first         : Load::load_v)       first->assemble();
//
//*********************************************      support
//
  for(auto const&  first         : Support::support_v) first->assemble();
//
//*********************************************      
//
// lhs.dim(2,2);
// rhs.dim(2);
// for(int ii=0; ii<2; ii++) lhs(ii,ii) = 2.0;
// lhs(1,0) = 0;
// lhs(0,1) = 1;
// for(int ii=0; ii<2; ii++) rhs(ii)    = 1.0;
}     


//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Solver::solve  -----
//
//
// C: Lapack solver
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Solver::solve(MDouble& mat, ADouble& vec)
{
  bool is_transpose = true;       // transpose for c to fortran compatibility
  bool is_symmetric = true;       // info on upper triangle is required. 
                                  // if lower triangle is given, transpose is not necessary
  
  char uplo='U';
  int  info;
  int  nrhs  = 1;
  int  neq   = mat.size1();
  int  lda   = mat.maxsize2();
  int  ldb   = vec.maxsize();
  int    nrow = mat.size1();
  int    ncol = mat.size2();

  AInt ipivot; 
  ipivot.dim(neq);
//
//*********************************************   
//
  std::cout<<"        : rhs  :        "<<vec<<std::endl;
  std::cout<<"        : lhs  :        "<<std::endl;
  std::cout<<mat;
  
//
//*********************************************      self trasnpose
//
  if(is_transpose){
    double aux;

    if(nrow != ncol)
      diag_mesg(diag.error,"SlvDenseLU::solve: Matrix is not square:"<<std::endl);

    for (int ii = 0; ii < nrow; ii++)
      for (int jj = ii; jj < ncol; jj++)
      {
        aux = mat(ii, jj);
        mat(ii, jj) = mat(jj, ii);
        mat(jj, ii) = aux;
      }
  }

//
//*********************************************      symmetrix
//
  if(is_symmetric){
    dposv_(&uplo,&neq,&nrhs,mat.pointer(),&lda,               vec.pointer(),&ldb,&info); //symmetric. works
  }
//
//*********************************************      nonsymetrix
//
  else{
    dgesv_(     &neq,&nrhs,mat.pointer(),&lda,ipivot.pointer(),vec.pointer(),&ldb,&info);
  }
//
//*********************************************      
//
  {
    std::cout<<"        : guess:        "<<vec<<std::endl;
  }

  if(info > 0){
    diag_mesg(diag.error,"SlvDenseLu::dposv: Solver error. Error num:"<<info<<std::endl);
  }
}
//
//*********************************************      new from here. Real.
//

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Solver::add_rhs  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Solver::add_rhs(AInt& index, double val)
{
  for(int ii=0; ii<index.size(); ii++)
    rhs(index(ii)) += val;
}

void Solver::add_lhs(AInt& index, double val)
{
  for(int ii=0; ii<index.size(); ii++)
    lhs(index(ii),index(ii)) += val;
}
