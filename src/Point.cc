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

