#include "Point.hh"
#include "Load.hh"
#include "Plate.hh"
#include "Essential.hh"

void Plate::init()
{
  ninteg   = 2;
  area     = d_zero;
}

void Plate::SamplePoint(int integ){

}


void Plate::add_edge(){
}

void Plate::assemble(){
}

void Plate::potential(){

  for(int ii=0; ii<ninteg; ii++){

    SamplePoint(ii);
    Grad(ii);
    Force();
    Fint();
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

}

void Plate::ShapeH()
{
    
}

void Plate::Force(){}
void Plate::Fint(){}
void Plate::Stretch(){}
void Plate::Stiffness(){}
