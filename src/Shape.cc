#include "Plate.hh"
#include "Diagnostic.hh"

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Plate::integration  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Plate::integration_point()
{

  switch(ninteg){
    
  case 3:
    xg_v[0] = { 1./6., 2./3., 1./6. };
    xg_v[1] = { 1./6., 1./6., 2./3. };
    wg_v    = { 1./3. , 1./3., 1./3. };
    break;
    
  case 7:
    xg_v[0] = {0.1012865073235, 0.7974269853531, 0.1012865073235, 0.4701420641051, 0.4701420641051, 0.0597158717898, 0.3333333333333};
    xg_v[1] = {xg_v[0][0],      xg_v[0][0],      xg_v[0][1],      xg_v[0][5],      xg_v[0][3],      xg_v[0][3],      xg_v[0][6]     };
    wg_v    = {0.1259391805448, 0.1259391805448, 0.1259391805448, 0.1323941527885, 0.1323941527885, 0.1323941527885, 0.225};
    break;
    
  case 13:
    xg_v[0] = { 0.0651301029022, 0.8697397941956, 0.0651301029022, 0.3128654960049, 0.6384441885698, 0.0486903154253, 0.6384441885698, 0.3128654960049, 0.0486903154253, 0.2603459660790, 0.4793080678419, 0.2603459660790, 0.3333333333333 };
    xg_v[1] = { xg_v[0][0],      xg_v[0][0],      xg_v[0][1],      xg_v[0][5],      xg_v[0][3],      xg_v[0][4],      xg_v[0][5],      xg_v[0][4],      xg_v[0][3],      xg_v[0][9],      xg_v[0][9],      xg_v[0][10],     xg_v[0][12]};
    wg_v    = { 0.0533472356088, 0.0533472356088, 0.0533472356088, 0.0771137608903, 0.0771137608903, 0.0771137608903, 0.0771137608903, 0.0771137608903, 0.0771137608903, 0.1756152574332, 0.1756152574332, 0.1756152574332, -0.1495700444677};
    
    break;
  default:
    diag_mesg(diag.error,"Plate::integration_point. Using an ilegal number of integration points:"<<ninteg<<std::endl);
  }
  
  if(shape_origin == edge_t)
    for(auto& vv : xg_v[0]) vv -= 0.5;
  
  for(auto& ww : wg_v) ww *= 0.5;

  double ans;

  diag<<"ninteg:"<<ninteg<<std::endl;
  diag<<"shape_origin:"<<shape_origin<<std::endl;

  ans = d_zero;  
  for(auto vv : xg_v[0]){
    ans += vv;
  }
  diag<<"xg0:"<<ans<<std::endl;

  ans = d_zero;  
  for(auto vv : xg_v[1]){
    ans += vv;
  }
  diag<<"xg1:"<<ans<<std::endl;

  ans = d_zero;    
  for(auto vv : wg_v){
    ans += vv;
  }
  diag<<"wg:"<<ans<<std::endl;

}


//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Plate::Shape  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Plate::Shape()
{
  switch(shape_origin){
  case edge_t:
    ShapeEdge();
    break;
  case  corner_t:
    ShapeCorner();
    break;
  }
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  void Plate::ShapeCorner  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Plate::ShapeCorner()
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
//               -----  void Plate::ShapeEdge  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Plate::ShapeEdge()
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
