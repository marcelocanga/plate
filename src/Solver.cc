#include "Plate.hh"
#include "Load.hh"
#include "Solver.hh"

void Solver::assemble()
{
  //  for(auto const& [first,second] : Plate::element_m) second->build_index();
  for(auto const& [first,second] : Plate::plate_m) second->assemble();
  for(auto const&  first         : Load::load_v)       first->assemble();
}     

void Solver::solve()
{}
