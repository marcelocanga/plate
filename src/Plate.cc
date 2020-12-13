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
MDouble Plate::b_shape, Plate::b_shape_h;
  
ADouble Plate::fint;
MDouble Plate::stiff;

void Plate::init()
{
  ninteg   = 3;
  area     = d_zero;

  nnode    = 3;
  nedge    = 3;

  nidof    = 2;            
  eldim    = 2;

  ndof     = nnode*eldim +nedge+nidof;

  shape.     resize(4);
  shape_h.   resize(3);
  b_shape.   resize(4,2);
  b_shape_h. resize(3,2);

  edof_loc.  resize(ndof);
  fint.      resize(ndof);
  stiff.     resize(ndof,ndof);
}

void Plate::count_dof(int& icount)
{

  for(int jj=0;jj<nnode;jj++)
      for(int kk=0;kk<eldim;kk++)
	      edof_loc( jj*eldim + kk) = point_v[jj]->pcount_dof(kk,icount);

  for(int jj=0; jj<nedge; jj++)
      edof_loc(nnode*eldim+jj) = point_v[jj]->pcount_dof(0,icount);
  
  for(int jj=0; jj<nidof; jj++) edof_loc(nnode*eldim + nedge +jj) = ++icount;
}


void Plate::SamplePoint(int integ){
  wgt = wg [integ];
  xr  = xg[integ][0];
  xs  = xg[integ][1];
}


void Plate::add_edge(){

}

void Plate::assemble(){
}

void Plate::potential(){

  compute_constitutive();

  for(int ii=0; ii<ninteg; ii++){

    SamplePoint(ii);
    Grad(ii);
    Stiffness();
  }

}

void Plate::compute_constitutive()
{

}

void Plate::Grad(int integ)
{
    Shape();
    ShapeH();
    Stretch();
}

void Plate::Shape()
{
  shape(0) = 1.0 - xr - xs;
  shape(1) = xr;
  shape(2) = xs;
  shape(3) = xr*xs - xr*xr*xs - xs*xs*xr;
    
  b_shape(0,0) = -1.0; 
  b_shape(1,0) =  1.0;
  b_shape(2,0) =  0.0;
  b_shape(3,0) =  xs - 2 * xr * xs - xs * xs;
  
  b_shape(0,1) = -1.0;
  b_shape(1,1) =  0.0;
  b_shape(2,1) =  1.0;
  b_shape(3,1) =  xr -  xr * xr  - 2.0 * xr * xs;
}

void Plate::ShapeH()
{
  shape_h(0) =  1.0 - 2.0 * xs;
  shape_h(1) = -1.0 + 2.0 * ( xr + xs);
  shape_h(2) =  1.0 - 2.0 * xr;
    
  b_shape_h(0,0) =  0.0;
  b_shape_h(1,0) =  2.0 * xs; 
  b_shape_h(2,0) = -2.0;

  b_shape_h(0,1) = -2.0;
  b_shape_h(1,1) =  2.0 * xr;
  b_shape_h(2,1) =  0.0;

}

void Plate::Fint(){}
void Plate::Stretch(){}
void Plate::Stiffness(){
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
}
