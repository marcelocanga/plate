#include "Point.hh"
#include "Load.hh"
#include "Plate.hh"

void Plate::init()
{
  ninteg   = 2;
  tickness = 0.1;
  young    = 2e5;
  poisson  = 0.3;
}
void Plate::SamplePoint(int integ){

}


void Plate::add_point(){
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
