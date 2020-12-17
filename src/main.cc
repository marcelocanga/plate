#include<iostream>
#include<fstream>
#include<sstream>
#include "Point.hh"
#include "Load.hh"
#include "Plate.hh"
#include "Solver.hh"
#include "Support.hh"
#include "Report.hh"

std::map<std::string, Point*>    Point::point_m;
std::set<Point*,LtPoint>         Point::u_point_s;
std::map<std::string, Plate*>    Plate::plate_m;
std::vector<Load*>               Load::load_v;
std::vector<Support*>            Support::support_v;

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  int main  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

int main(int argc, char *argv[])
{
  std::string input_file;
  Solver sol;

  std::cout << std::fixed << std::showpoint;
  std::cout << std::setprecision(6);
  std::clog << std::fixed << std::showpoint;
  std::clog << std::setprecision(3);
//
//*********************************************      first argument, input file name
//
  if(argc > 1) input_file = argv[1];
  else         input_file = "plate.inp";
//
//*********************************************      
//
  rep.set_fos(input_file);
//
//*********************************************      parse input
//
  sol.parse_input(input_file);
//
//*********************************************      
//
  rep.fos_close();
//
//*********************************************      done
//
}
