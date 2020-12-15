#include<sstream>
#include "Plate.hh"
#include "Point.hh"

Point::Point(){
  dof_loc = -1;
}

std::string Point::new_name()
{
  std::stringstream ss;
  static int index;
  ss << index;
  index++;
  return "edge-"+ ss.str();
}
