#ifndef _DIAGNOSTIC_FEM_
#define _DIAGNOSTIC_FEM_

class Diagnostic;

extern Diagnostic diag;

class Diagnostic{

public:
  
  enum DiagLevel { none, error, info, echo, debug } diag_level;
  
  Diagnostic();
  bool level(enum DiagLevel);
  void set_level(int);
};


#endif

  
