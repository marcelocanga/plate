#include "Plate.hh"
#include "Load.hh"
#include "Solver.hh"
#include "Support.hh"
#include "Diagnostic.hh"

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
// C: Marcelo Canga. Dec 2020
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
         poisson_a, thickness_a, end_a, line_force_a, line_moment_a, diagnostic_a };
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
    else if (line.find("*lfor")   != std::string::npos) token = line_force_a;
    else if (line.find("*lmom")   != std::string::npos) token = line_moment_a;
    else if (line.find("*pres")   != std::string::npos) token = pressure_a;
    else if (line.find("*supp")   != std::string::npos) token = support_a;
    else if (line.find("*solv")   != std::string::npos) token = solve_a;    
    else if (line.find("*youn")   != std::string::npos) token = young_a;    
    else if (line.find("*pois")   != std::string::npos) token = poisson_a;    
    else if (line.find("*thic")   != std::string::npos) token = thickness_a;    
    else if (line.find("*diag")   != std::string::npos) token = diagnostic_a;
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
      diag_m(diag.info,"Point:"<<po_pt->name<<":"<<po_pt->coor(0)<<","<<po_pt->coor(1)<<","<<po_pt->coor(2)<<std::endl);
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
          diag_m(diag.error, "Point in plate not found:" << nam << ", plate" << el_pt->name << std::endl);
          break;
        }
      }
      Plate::plate_m[el_pt->name] = el_pt;

      el_pt->poisson   = def_poisson;
      el_pt->young     = def_young;
      el_pt->thickness = def_thickness;

      diag_m(diag.info,"Plate:"<<el_pt->name<<","<<point[0]<<","<<point[1]<<","<<point[2]<<std::endl);
    }
      break;
//.............................................      shear force
//
    case line_force_a:
    case force_a:{
      std::string ename;
      Load *lo_pt;
      if (token == force_a) lo_pt = new Load(Load::force_a);
      else                  lo_pt = new Load(Load::line_force_a);
      iss >> ename >>lo_pt->eside >> lo_pt->value;

      if(Plate::plate_m.count(ename))
        lo_pt->el_pt = Plate::plate_m[ename];
      else{
        diag_m(diag.error,"Plate in force not found:" << ename << std::endl);
        break;
      }

      Load::load_v.push_back(lo_pt);

      if (token == force_a){
        diag_m(diag.info,"Force. Plate:"<<ename <<", side:"<<lo_pt->eside <<", value:"<<lo_pt->value <<std::endl);
      }
      else{
        diag_m(diag.info,"Line Force. Plate:"<<ename <<", side:"<<lo_pt->eside <<", value:"<<lo_pt->value <<std::endl);
      }
    }
      break;
//.............................................      moment
//
    case moment_a:{
      int gdir;
      std::string pname;
      Load *lo_pt = new Load(Load::moment_a);
      
      iss >> pname >> lo_pt->gdir >> lo_pt->value;

      if(Point::point_m.count(pname))
        lo_pt->po_pt = Point::point_m[pname];
      else{
        diag_m(diag.error,"Point in moment not found:" << pname << std::endl);
        break;
      }

      Load::load_v.push_back(lo_pt);

      diag_m(diag.info,"Moment. Point:"<<pname <<", dir:"<<lo_pt->gdir<<", value:"<<lo_pt->value<<std::endl);
}
      break;
//.............................................      line moment
//
    case line_moment_a:{
      int gdir;
      std::string ename;
      Load *lo_pt = new Load(Load::line_moment_a);
      
      iss >> ename >>lo_pt->eside >> lo_pt->gdir >> lo_pt->value;

      if(Plate::plate_m.count(ename))
        lo_pt->el_pt = Plate::plate_m[ename];
      else{
        diag_m(diag.error,"Plate in line moment not found:" << ename << std::endl);
        break;
      }

      Load::load_v.push_back(lo_pt);

      diag_m(diag.info,"Line Moment. Plate:"<<ename <<", side:"<<lo_pt->eside <<", dir:"<<lo_pt->gdir<<", value:"<<lo_pt->value<<std::endl);
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
        diag_m(diag.error,"Plate in pressure not found:" << ename << std::endl);
        break;
      }
      diag_m(diag.info,"Plate Pressure in:"<<ename<<", val:"<<val<<std::endl);
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
        diag_m(diag.error,"Plate in support not found:" << ename << std::endl);
        break;
      }

      Support::support_v.push_back(su_pt);

      diag_m(diag.info,"Fix Support. Plate:"<<ename <<", side:"<<su_pt->eside <<std::endl);
    }
      break;
//.............................................      material properties/thickness
//
    case thickness_a:
      iss >> def_thickness;
      diag_m(diag.info,"Plate thickness:"<<def_thickness<<std::endl);
      break;
      
    case young_a:
      iss >> def_young;
      diag_m(diag.info,"Young's modulus:"<<def_young<<std::endl);
      break;

    case poisson_a:
      iss >> def_poisson;
      diag_m(diag.info,"Poisson's ratio:"<<def_poisson<<std::endl);
      break;
//.............................................      
//
    case diagnostic_a:{
      int level;
      iss >> level;
      diag.set_level(level);
    }
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
  diag_m(diag.info,"Solver starts."<<std::endl);
//
//*********************************************     add nodes at element edges 
//
  for(auto const& [first,second] : Plate::plate_m) second->add_edge();
//
//*********************************************      count and map equation index
//
  for(auto const& [first,second] : Plate::plate_m) second->count_dof(ndof);
  diag_m(diag.info,"Ndof:"<<ndof<<std::endl);
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
  for(auto const& [first,second] : Plate::plate_m) second->compute_stress();
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
//*********************************************      done
//
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
  bool is_symmetric = false;       // info on upper triangle is required. 
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
diag_l(diag.echo,
  diag<<"        : rhs  :        "<<vec<<std::endl;
  diag<<"        : lhs  :        "<<std::endl;
  diag<<mat;
);
//
//*********************************************      self trasnpose
//
  if(is_transpose){
    double aux;

    if(nrow != ncol)
      diag_m(diag.error,"SlvDenseLU::solve: Matrix is not square:"<<std::endl);

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
  diag_l(diag.echo,
    diag<<"        : guess:        "<<vec<<std::endl;
  );

  if(info > 0){
    diag_m(diag.error,"SlvDenseLu::dposv: Solver error. Error num:"<<info<<std::endl);
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

void Solver::set_ans(AInt& index, double val)
{
  for(int ii=0; ii<index.size(); ii++){
    for(int kk=0; kk<ndof; kk++){
      lhs(kk,index(ii)) = d_zero;
      lhs(index(ii),kk) = d_zero;
    }
    lhs(index(ii),index(ii)) = d_one;
    rhs(index(ii))           = val;
  }
}
