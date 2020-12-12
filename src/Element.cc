#include "Plate.hh"

void Element::init()
{
  ninteg   = 2;
  tickness = 0.1;
  young    = 2e5;
  poisson  = 0.3;
}
void Element::SamplePoint(int integ){

}


void Element::add_point(){
}

void Element::assemble(){
}

void Element::potential(){

  for(int ii=0; ii<ninteg; ii++){

    SamplePoint(ii);
    Grad(ii);
    Force();
    Fint();
    Stiffness();
  }

}

void Element::compute_constitutive()
{

}

void Element::Grad(int integ)
{
    Shape();
    ShapeH();
    Stretch();
}

void Element::Shape()
{

}

void Element::ShapeH()
{
    
}

void Element::Force(){}
void Element::Fint(){}
void Element::Stretch(){}
void Element::Stiffness(){}