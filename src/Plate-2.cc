#include "Point.hh"
#include "Load.hh"
#include "Plate.hh"
#include "Essential.hh"
#include "Array.hh"
#include "Solver.hh"
#include "Diagnostic.hh"

double Plate::wg[3]={1.0/6.0,1.0/6.0,1.0/6.0};
double Plate::xg[2][3]={-1./3.,1./3.,0.,
                        1./6.,1./6.,2./3.};
//double Plate::xg[2][3]={-1./2.,1./2., 0.0,
//                          0.,  0.0,   1.0};
//double Plate::xg[2][3]={  0., 1./4., -1./4.,
//                         0., 1./2.,  1./2.};

int Plate::ip;
double Plate::wgt,Plate::xr,Plate::xs;
ADouble Plate::shape,Plate::shape_h;
MDouble Plate::d_shape, Plate::d_shape_h;
MDouble Plate::b_grad, Plate::w_grad;
MDouble Plate::constitutive_b;
MDouble Plate::constitutive_s;
ADouble Plate::fint;
MDouble Plate::stiff;

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//
//               -----  Plate::Plate  -----
//
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Plate::Plate(){
  init();
  area     = d_zero;
  d_area   = d_zero;
  pressure = d_zero;
}


//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Plate::init  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Plate::init()
{
  ninteg   = 3;
  nnode    = 3;
  nedge    = 3;
  
  nidof    = 1;
  nshear   = 1;
  
  eldim    = 2;

  nedof     = (nnode+nidof+nshear)*eldim + nedge ;

  shape.     dim(4);
  shape_h.   dim(3);
  d_shape.   dim(4,2);
  d_shape_h. dim(3,2);
  b_grad.    dim(3,8);
  w_grad.    dim(2,3);

  fint.      dim(nedof);
  stiff.     dim(nedof,nedof);

  constitutive_b.dim(3,3);
  constitutive_s.dim(2,2);
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Plate::get_index  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

bool Plate::get_index(int side, int idir, enum IndexType type, AInt& index)
{
  int nd1     =   side;
  int nd2     =  (side+1)%3;

  switch(type){
  case support_t:
    index.resize(5);
    index(0) = edof_loc[2*nd1];
    index(1) = edof_loc[2*nd1+1];
    index(2) = edof_loc[2*nd2];
    index(3) = edof_loc[2*nd2+1];
    index(4) = edof_loc[(nnode+nidof+nshear)*eldim+nd1];
    diag_m(diag.echo,"Plate::get_index:support:"<<index<<std::endl);
    break;
  case moment_t:
    index.resize(2);
    index(0) = edof_loc[2*nd1+idir];
    index(1) = edof_loc[2*nd2+idir];
    diag_m(diag.echo,"Plate::get_index:moment:"<<index<<std::endl);
    break;
  case force_t:
    index.resize(1);
    index(0) = edof_loc[(nnode+nidof+nshear)*eldim+nd1];
    diag_m(diag.echo,"Plate::get_index:force:"<<index<<std::endl);
    break;
  }

  return true;
  
}


//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Plate::count_dof  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Plate::count_dof(int& ndof)
{
  Point* po_pt;

  for(int jj=0;jj<nnode;jj++){
    po_pt = point_v[jj];
    if(po_pt->dof_loc < 0){ 
      po_pt->dof_loc = ndof;
      ndof +=2;
    }
    edof_loc.push_back(po_pt->dof_loc);
    edof_loc.push_back(po_pt->dof_loc+1);
  }
  
  for(int ii=0; ii<(nshear+nidof)*eldim; ii++)
    edof_loc.push_back(ndof+ii);
  ndof += (nshear+nidof)*eldim;

  for(int jj=0;jj<nedge;jj++){
    po_pt = point_v[jj+nnode];
    if(po_pt->dof_loc<0){
      po_pt->dof_loc = ndof;
      ndof += 1;
    }
    edof_loc.push_back(po_pt->dof_loc);
  }

diag_l(diag.echo,
std::clog<<"Element dof, name:"<<name<<std::endl;
for(int ii=0; ii<13; ii++)
std::clog<<edof_loc[ii]<<" ";
std::clog<<std::endl;
);

}




//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Plate::add_edge  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Plate::add_edge(){
  
  Coord coor;
  Point  po_aux;
  Point *po_pt;

  std::clog<<"Add Edge: Plate:"<<name<<std::endl;

  for(int ii=0; ii<nnode; ii++){                     // check if there is point with this coord
    int iie = (ii+1) % nnode;
    scale(0.5,point_v[ii]->coor,        po_aux.coor);
    scale(0.5,point_v[iie]->coor, d_one,po_aux.coor);

    if(Point::u_point_s.count(&po_aux)){ 
      auto it = Point::u_point_s.find(&po_aux);
      point_v.push_back(*it);
      diag_m(diag.echo,"Already found edge: side:"<<iie<<", name:"<<(*it)->name<<std::endl);
    }
    else{
      po_pt = new Point();
      po_pt->name = Point::new_name();
      scale(po_aux.coor,po_pt->coor);
      point_v.push_back(po_pt);
      Point::u_point_s.insert(po_pt);
      diag_m(diag.echo,"New edge: side"<<iie<<", name:"<<po_pt->name<<","<<po_pt->coor<<std::endl);
    }

  }

}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Plate::assemble  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Plate::assemble(){

  potential();

  for(int ii=0; ii<nedof; ii++){
    for(int jj=0; jj<nedof; jj++)
      Solver::current()->lhs(edof_loc[ii],edof_loc[jj]) += stiff(ii,jj);
    Solver::current()->rhs(edof_loc[ii])                += fint(ii);
  }
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Plate::potential  -----
//
//
// C: build stiffness matrix
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Plate::potential(){
  area  = d_zero;
  fint  = d_zero;
  stiff = d_zero;

  compute_constitutive();

  for(int ii=0; ii<ninteg; ii++){

    SamplePoint(ii);
    Grad(ii);
    Fint();
    Stiffness();
  }

diag_l(diag.echo,
  std::clog<<"Area:"<<area<<std::endl;
  std::clog<<"Fint:"<<fint<<std::endl;
);
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Plate::stress  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Plate::compute_stress(){
  area  = d_zero;

  compute_constitutive();

  for(int ii=0; ii<ninteg; ii++){

    SamplePoint(ii);
    Grad(ii);
    Stress();
  }

diag_l(diag.echo,
  std::clog<<"Area:"<<area<<std::endl;
  std::clog<<"Fint:"<<fint<<std::endl;
);
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Plate::SamplePoint  -----
//
//
// C: isoparametric coords and weights
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Plate::SamplePoint(int integ){
  wgt = wg [integ];
  xr  = xg[0][integ];
  xs  = xg[1][integ];

diag_l(diag.debug,
  std::clog<<"integ:"<<integ<<std::endl;
  std::clog<<"wgt:"<<wgt<<", xr:"<<xr<<", xs:"<<xs<<std::endl;
);
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Plate::compute_constitutive  -----
//
//
// C: bending and shear elasticity matrix
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Plate::compute_constitutive()
{
  double kfac = d_one;
  double fac1 = young*thickness*thickness/(d_one-poisson*poisson);
  double fac2 = young*kfac/ 2.0 / (d_one+poisson);

  constitutive_b = d_zero;

  constitutive_b(0,0) = fac1;
  constitutive_b(0,1) = fac1 * poisson; 
  constitutive_b(1,0) = fac1 * poisson;
  constitutive_b(1,1) = fac1;
  constitutive_b(2,2) = fac1 * (d_one - poisson) / 2.0;

  constitutive_s      = d_zero;
  constitutive_s(0,0) = fac2;
  constitutive_s(1,1) = fac2;

diag_l(diag.debug,
  std::clog<<"Constitutive: b"<<constitutive_b<<std::endl;
  std::clog<<"Constitutive: s"<<constitutive_s<<std::endl;
);

}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Plate::Shape  -----
//
//
// C: conforming and non-conforming triangle base functions
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Plate::Shape()
{
  shape(0) = 0.5 * (1.0 - 2.0 * xr - xs );
  shape(1) = 0.5 * (1.0 + 2.0 * xr - xs );
  shape(2) = xs;
  shape(3) = 0.25*xs - xr*xr*xs - 0.5*xs*xs + 0.25*xs*xs*xs;
    
  d_shape(0,0) = -1.0; 
  d_shape(1,0) =  1.0;
  d_shape(2,0) =  0.0;
  d_shape(3,0) =  -2.0 * xr * xs;
  
  d_shape(0,1) = -0.5;
  d_shape(1,1) = -0.5;
  d_shape(2,1) =  1.0;
  d_shape(3,1) =  0.25 - xr* xr - xs + 0.75 * xs * xs;

  shape_h(0) =  1.0       - 2.0 * xs;
  shape_h(1) =  2.0 * xr  +       xs;
  shape_h(2) = -2.0 * xr  +       xs;
    
  d_shape_h(0,0) =  0.0;
  d_shape_h(1,0) =  2.0; 
  d_shape_h(2,0) = -2.0;

  d_shape_h(0,1) = -2.0;
  d_shape_h(1,1) =  1.0;
  d_shape_h(2,1) =  1.0;
}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Plate::Grad  -----
//
//
// C: gradients
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Plate::Grad(int integ)
{
  M2Double dx, dxi;
  
  Shape();

//
//*********************************************     Jacobian, Dx/Dr, d_area, Dr/Dx
//
   for(int ii=0; ii<2; ii++)
    for(int jj=0; jj<2; jj++){
      double fac = d_zero;
      for(int kk=0; kk<nnode; kk++) fac += d_shape(kk,jj) * point_v[kk]->coor(ii); 
      dx(ii,jj) = fac;
    }

  d_area = dx(0,0) * dx(1,1) - dx(1,0) * dx(0,1);

  dxi(1,1) =  dx(0,0) / d_area;
  dxi(0,0) =  dx(1,1) / d_area;
  dxi(0,1) = -dx(0,1) / d_area;
  dxi(1,0) = -dx(1,0) / d_area;
//
//*********************************************      bending, Exx,Eyy,Exy
//
  b_grad = d_zero;
  
  for(int kk=0; kk<4; kk++)
    for(int jj =0; jj<2; jj++){
      b_grad(0,2*kk+0) +=       d_shape(kk,jj) * dxi(jj,0) ;
      b_grad(1,2*kk+1) +=       d_shape(kk,jj) * dxi(jj,1) ;
      b_grad(2,2*kk+1) += 0.5 * d_shape(kk,jj) * dxi(jj,0);      
      b_grad(2,2*kk+0) += 0.5 * d_shape(kk,jj) * dxi(jj,1);      
    }

//
//*********************************************      Wx,Wy
//
  w_grad = d_zero;

  for(int kk=0; kk<3; kk++)
    for(int jj =0; jj<2; jj++){
      w_grad(0,kk) +=  d_shape_h(kk,jj) * dxi(jj,0);
      w_grad(1,kk) +=  d_shape_h(kk,jj) * dxi(jj,1);      
    }
//
//*********************************************      
//
diag_l(diag.debug,
  std::clog<<"d_area:"<<d_area<<std::endl;
  std::clog<<"dx"<<std::endl;
  std::clog<<dx<<std::endl;

  std::clog<<"dxi"<<std::endl;
  std::clog<<dxi<<std::endl;

  std::clog<<"Shape"<<std::endl;
  std::clog<<shape<<std::endl;

  std::clog<<"Shape-h"<<std::endl;
  std::clog<<shape_h<<std::endl;

  std::clog<<"d-Shape"<<std::endl;
  std::clog<<d_shape<<std::endl;

  std::clog<<"d-Shape-h"<<std::endl;
  std::clog<<d_shape_h<<std::endl;


  std::clog<<"Curvature"<<std::endl;
  std::clog<<b_grad<<std::endl;

  std::clog<<"w-gradient"<<std::endl;
  std::clog<<w_grad<<std::endl;
);
//
//*********************************************      done
//
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Plate::Fint  -----
//
//
// C: pressure forces on the plate face, and area calc
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Plate::Fint(){
  int n2     =  eldim*(nnode+nidof+nshear);
  double fac =  d_area * wgt ;    
  
  area          += fac ;
  
  for(int mm=0; mm<nedge; mm++){
    fint(n2 + mm) += fac * pressure * shape_h(mm) ;
  }
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Plate::Stiffness  -----
//
//
// C: stiffness matrix, 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Plate::Stiffness(){

  double fac = d_area * wgt * thickness;

//
//*********************************************      bending 
//
  int ncons = nnode+nidof;
  
  for(int nn=0; nn<ncons*eldim; nn++)
    for(int mm=0; mm<ncons*eldim; mm++)
      for(int ii=0; ii<3; ii++)
	for(int jj=0; jj<3; jj++)
	  stiff(nn,mm)                          += fac * b_grad(ii,nn)*constitutive_b(ii,jj)*b_grad(jj,mm) ;

//
//*********************************************      shear
//
  for(int nn=0; nn<2; nn++)
    stiff(ncons*eldim+nn,ncons*eldim+nn)        += fac * constitutive_s(nn,nn);
//
//*********************************************      constraints shear,bending
//
  for(int mm=0; mm<ncons; mm++)
      for(int kk=0; kk<nshear*eldim; kk++){
	stiff(ncons*eldim+kk,mm*2+kk)            += fac * shape(mm);
	stiff(mm*2+kk,           ncons*eldim+kk) += fac * shape(mm);
      }
//
//*********************************************      constrains shear, grad displ
//

  int n3 = eldim*ncons;
  int n2 = eldim*(ncons+nshear);
  
  for(int mm=0; mm<nedge; mm++)
    for(int kk=0; kk<2; kk++){
      stiff(n3+kk,mm+n2)         -= fac * w_grad(kk,mm);
      stiff(mm+n2,        n3+kk) -= fac * w_grad(kk,mm);
    }
//
//*********************************************      done
//
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Plate::Stress  -----
//
//
// C: Compute stresses at integration points
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Plate::Stress(){
  
  ADouble bending(8);
  ADouble shear(2);

  int ncons = nnode+nidof;
//
//*********************************************      
//
  for(int ii=0; ii<ncons*eldim; ii++)
    bending(ii) = Solver::current()->rhs(ii);
  for(int ii=0; ii<2; ii++)
    bending(ii) = Solver::current()->rhs(ncons*eldim+ii);  
//
//*********************************************      bending 
//
  for(int mm=0; mm<ncons*eldim; mm++)
    for(int ii=0; ii<3; ii++)
      for(int jj=0; jj<3; jj++)
	bending(ii)   += thickness * constitutive_b(ii,jj)*b_grad(jj,mm)*bending(mm) ;
//
//*********************************************      shear
//
  for(int nn=0; nn<2; nn++)
    shear(nn)         += thickness *  constitutive_s(nn,nn)*shear(nn);
//
//*********************************************      done
//
}

  
//    for(int mm=0; mm<nnode*eldim; mm++)
//      for(int kk=0; kk<nalfa; kk++){
//	stiff(nnode*eldim+kk,mm)            += shapeh(kk)*d_stretch[0](mm) *iso_jacob*wgt;
//	stiff(mm,           nnode*eldim+kk) += shapeh(kk)*d_stretch[0](mm) *iso_jacob*wgt;
//      }


//      for(int mm=0; mm<nalfa; mm++)
//	for(int kk=0; kk<nalfa; kk++)
//	  stiff(nnode*eldim+kk,nnode*eldim+mm) -= shapeh(kk)*shapeh(mm)*elastic*iso_jacob*wgt;

  
//    for(int mm=0; mm<nnode*eldim; mm++)
//      for(int kk=0; kk<nalfa; kk++){
//	stiff(nnode*eldim+kk,mm)            += shapeh(kk)*d_stretch[0](mm) *iso_jacob*wgt;
//	stiff(mm,           nnode*eldim+kk) += shapeh(kk)*d_stretch[0](mm) *iso_jacob*wgt;
//      }

//    double elastic = d_one/constitutive(0,0);
// 
//  for(int nn=0; nn<nnode*eldim; nn++)
//    for(int mm=0; mm<nnode*eldim; mm++)
//      for(int ii=1; ii<eldim; ii++)
//	for(int jj=1; jj<eldim; jj++)
//	  stiff(nn,mm)                      += d_stretch[ii](nn)*constitutive(ii,jj)*d_stretch[jj](mm) *iso_jacob*wgt;
//
//    for(int mm=0; mm<nnode*eldim; mm++)
//      for(int pp=0; pp<nnode*eldim; pp++)
//	for(int ii=0; ii<eldim; ii++)
//	  stiff(mm,pp)                      += force(ii)*dd_stretch[ii](mm,pp)*iso_jacob*wgt;
//
//    for(int mm=0; mm<nnode*eldim; mm++)
//      for(int kk=0; kk<nalfa; kk++){
//	stiff(nnode*eldim+kk,mm)            += shapeh(kk)*d_stretch[0](mm) *iso_jacob*wgt;
//	stiff(mm,           nnode*eldim+kk) += shapeh(kk)*d_stretch[0](mm) *iso_jacob*wgt;
//      }
//
//    if(is_hybrid_elastic)
//      for(int mm=0; mm<nalfa; mm++)
//	for(int kk=0; kk<nalfa; kk++)
//	  stiff(nnode*eldim+kk,nnode*eldim+mm) -= shapeh(kk)*shapeh(mm)*elastic*iso_jacob*wgt;
//

