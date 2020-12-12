#include "Plate.hh"


void Solution::assemble()
{
  //  for(auto const& [first,second] : Element::element_m) second->build_index();
  for(auto const& [first,second] : Element::element_m) second->assemble();
  for(auto const&  first         : Load::load_v)       first->assemble();
}     

void Solution::solve()
{}
