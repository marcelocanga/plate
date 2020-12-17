#ifndef _REPORT_FEM_
#define _REPORT_FEM_

class Report;

extern Report rep;

class Report{

  std::ofstream fos;
  
public:
  
  Report();
  void set_fos(std::string);
  void fos_close();
  void summary();
//
//************************************************      operators
//
  template< class T>

  Report& operator<<(T& p){
    fos<<p;
    std::clog<<p;
    return *this; 
  }
  Report& operator<<(std::ostream& omanip(std::ostream &));

};


#endif

  
