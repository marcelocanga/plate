#ifndef _DIAGNOSTIC_FEM_
#define _DIAGNOSTIC_FEM_

class Diagnostic;

extern Diagnostic diag;

class Diagnostic{

public:
  
  enum DiagLevel { none, error, info, echo, debug, detail } diag_level;
  
  Diagnostic();
  bool level(enum DiagLevel);
  void set_level(int);
  std::string top_msg(enum DiagLevel level);
//
//************************************************      operators
//
  template< class T>

  Diagnostic& operator<<(T& p){
    std::clog<<p;
    return *this; 
  }
  Diagnostic& operator<<(std::ostream& omanip(std::ostream &));

};


#endif

  
