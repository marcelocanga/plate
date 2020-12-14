#include "Point.hh"
#include "Load.hh"
#include "Plate.hh"
#include "Essential.hh"
#include "Array.hh"


double Plate::wg[3]={1.0/3.0,1.0/3.0,1.0/3.0};
double Plate::xg[3][2]={1./6.,2./3.,1./6.,
                        1./6.,1./6.,2./3.};
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
  static bool is_first = true;

  if(is_first){
    init();
    is_first = false;
  }

  area   = d_zero;
  d_area = d_zero;

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

  ndof     = (nnode+nidof+nshear)*eldim + nedge ;

  shape.     dim(4);
  shape_h.   dim(3);
  d_shape.   dim(4,2);
  d_shape_h. dim(3,2);
  b_grad.    dim(3,8);
  w_grad.    dim(2,3);

  edof_loc.  dim(ndof);
  fint.      dim(ndof);
  stiff.     dim(ndof,ndof);

  constitutive_b.dim(3,3);
  constitutive_s.dim(2,2);
}


//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Plate::count_dof  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Plate::count_dof(int& icount)
{

  for(int jj=0;jj<nnode;jj++)
      for(int kk=0;kk<eldim;kk++)
	      edof_loc( jj*eldim + kk) = point_v[jj]->pcount_dof(kk,icount);

  for(int jj=0; jj<nedge; jj++)
      edof_loc(nnode*eldim+jj) = point_v[jj]->pcount_dof(0,icount);
  
  for(int jj=0; jj<nidof; jj++) edof_loc(nnode*eldim + nedge +jj) = ++icount;
}




void Plate::add_edge(){

}

void Plate::assemble(){
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

  compute_constitutive();

  for(int ii=0; ii<ninteg; ii++){

    SamplePoint(ii);
    Grad(ii);
    Stiffness();
  }

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
  xr  = xg[integ][0];
  xs  = xg[integ][1];
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
  constitutive_b(0,1) = fac1 * poisson ; 
  constitutive_b(1,0) = fac1;
  constitutive_b(1,1) = fac1 * poisson;
  constitutive_b(2,2) = fac1 * (d_one - poisson) / 2.0;

  constitutive_s      = d_zero;
  constitutive_s(0,0) = fac2;
  constitutive_s(1,1) = fac2;

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
  shape(0) = 1.0 - xr - xs;
  shape(1) = xr;
  shape(2) = xs;
  shape(3) = xr*xs - xr*xr*xs - xs*xs*xr;
    
  d_shape(0,0) = -1.0; 
  d_shape(1,0) =  1.0;
  d_shape(2,0) =  0.0;
  d_shape(3,0) =  xs - 2.0 * xr * xs - xs * xs;
  
  d_shape(0,1) = -1.0;
  d_shape(1,1) =  0.0;
  d_shape(2,1) =  1.0;
  d_shape(3,1) =  xr - 2.0 * xr * xs - xr * xr ;

  shape_h(0) =  1.0 - 2.0 * xs;
  shape_h(1) = -1.0 + 2.0 * ( xr + xs);
  shape_h(2) =  1.0 - 2.0 * xr;
    
  d_shape_h(0,0) =  0.0;
  d_shape_h(1,0) =  2.0; 
  d_shape_h(2,0) = -2.0;

  d_shape_h(0,1) = -2.0;
  d_shape_h(1,1) =  2.0;
  d_shape_h(2,1) =  0.0;
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
  M2Double xj, xji;
  
  Shape();

//
//*********************************************     Jacobian, Dx/Dr, d_area, Dr/Dx
//
   for(int ii=0; ii<2; ii++)
    for(int jj=0; jj<2; jj++){
      double fac = d_zero;
      for(int kk=0; kk<nnode; kk++) fac += d_shape(kk,jj) * point_v[kk]->coor(ii); 
      xj(ii,jj) = fac;
    }

  d_area = xj(0,0) * xj(1,1) - xj(1,0) * xj(0,1);

  xji(1,1) = xj(0,0)/ d_area;
  xji(0,0) = xj(1,1)/ d_area;
  xji(0,1) =-xj(0,1)/ d_area;
  xji(1,1) =-xj(1,0)/ d_area;
//
//*********************************************      bending, Exx,Eyy,Exy
//
  b_grad = d_zero;
  
  for(int kk=0; kk<4; kk++)
    for(int jj =0; jj<2; jj++){
      b_grad(0,2*kk+0) += d_shape(kk,jj) * xji(jj,0) ;
      b_grad(1,2*kk+1) += d_shape(kk,jj) * xji(jj,1) ;
      b_grad(2,2*kk+0) += 0.5 * d_shape(kk,jj) * xji(jj,0);      
      b_grad(2,2*kk+1) += 0.5 * d_shape(kk,jj) * xji(jj,1);      
    }

//
//*********************************************      Wx,Wy
//
  w_grad = d_zero;

  for(int kk=0; kk<3; kk++)
    for(int jj =0; jj<2; jj++){
      w_grad(0,kk+0) +=  d_shape_h(kk,jj) * xji(jj,0);
      w_grad(1,kk+1) +=  d_shape_h(kk,jj) * xji(jj,1);      
    }
//
//*********************************************      done
//
}

void Plate::Fint(){}
void Plate::Stretch(){}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Plate::Stiffness  -----
//
//
// C: stiffness matrix
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
	stiff(ncons*eldim+kk,mm*2+0)            += fac * shape(mm);
	stiff(ncons*eldim+kk,mm*2+1)            += fac * shape(mm);
	stiff(mm*2+0,           ncons*eldim+kk) += fac * shape(mm);
	stiff(mm*2+1,           ncons*eldim+kk) += fac * shape(mm);
      }
//
//*********************************************      constrains shear, grad displ
//

  int n2 = eldim*(ncons+nshear);
  
  for(int mm=0; mm<nedge; mm++)
    for(int kk=0; kk<2; kk++){
      stiff(n2+kk,mm)            -= fac * w_grad(kk,mm);
      stiff(mm,           n2+kk) -= fac * w_grad(kk,mm);
    }
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

