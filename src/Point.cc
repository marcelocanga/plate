#include<sstream>
#include "Plate.hh"
#include "Point.hh"


std::string Point::new_name()
{
  std::stringstream ss;
  static int index;
  ss << index;
  index++;
  return "edge-"+ ss.str();
}

int& Point::pcount_dof(const int& k,int& icount)
{
  if(pdof_loc(k)== -1){
    pdof_loc(k)=++icount;
    return icount;
  }
  else return pdof_loc(k);
}


void Point::pclean_dof()
{
  pdof_loc    = -1;
}
