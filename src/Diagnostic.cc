#include <iostream>
#include "Diagnostic.hh"
#include "Essential.hh"

Diagnostic diag;

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//
//               -----  Diagnostic::Diagnostic  -----
//
//
//
// C: Marcelo Canga. Dec 2020
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Diagnostic::Diagnostic()
{
  diag_level = echo;
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  bool Diagnostic::level  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

bool Diagnostic::level(enum DiagLevel level)
{
  if(level <= diag_level) return true;
  else                    return false;
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  bool Diagnostic::set_level  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void Diagnostic::set_level(int level)
{
  if(level < 0 || level > debug){
    diag_m(diag.error,"Diagnostic::set_level. Invalid level"<<std::endl);
    return;
  }
  diag_level = DiagLevel(level);
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  Diagnostic& Diagnostic::operator<<(std::ostream& omanip  -----
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Diagnostic& Diagnostic::operator<<(std::ostream& omanip(std::ostream &)){
      omanip(std::clog);
  return *this; 
}
