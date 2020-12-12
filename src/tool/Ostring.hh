#ifndef __OSTRING_MEC_HH__
#define __OSTRING_MEC_HH__

#include <iostream>
#include <sstream>
#include <iomanip>
#include <complex>
#include "String.hh"
#include "Map.hh"
#include "Vector.hh"
#include "Array.hh"


class Ostring : public std::ostringstream
{
public:
  enum Format {scientific,fixed,right,left,quote,none};
  Ostring();
  Ostring(double width,enum Format ft=none );

  void set_format(double width,enum Format ft); //10.2
  void set_format(double width); //10.2
  void set_format(enum Format ft);
  
  void set_format_plot(double value);
  void set_format_cpp_scientific();       // not used for double, produces standard c++, incosistent in windows

  enum Format get_format(){return format;}

  Ostring& put(bool);
  Ostring& put(int);
  Ostring& put(long int);
  Ostring& put(float);
  Ostring& put(double);
  Ostring& put(std::complex<double>);
  Ostring& put(std::string);
  Ostring& put(const char*);

  void clear(){ str(s_zero); }

//
//*********************************************    for formatted output:  copy from analyzer, notice the >> operator
//
  template <class T>
  Ostring& operator >> (T& p);
  template <class T>
  Ostring& operator >> (Vector<T>& p);
  template <class T>
  Ostring& operator >> (Array<T,1>& p);
  template <class T>
  Ostring& operator >> (Array<T,2>& p);
  template <class T, int Dim>
  Ostring& operator >> (Tiny<T,Dim>& p);

  Ostring& operator >> (const char* p);  // specialization to output cons
//
//*********************************************      
//

private:

  int width_;
  int precision_;
  enum Format format;

};

inline 
Ostring&  Ostring::operator >> (const char* p){
  return put(p);
}
template <class T> Ostring&  Ostring::operator >> (T& p){
  return put(p);
}
template <class T> Ostring& Ostring::operator >> (Vector<T>& p){
  for(int ii=0; ii < p.size(); ii++)
    put(p[ii]);
  return *this;
}
template <class T> Ostring& Ostring::operator >> (Array<T,1>& p){
  for(int ii=0; ii < p.size(); ii++)
    put(p(ii));
  return *this;
}

template <class T> Ostring& Ostring::operator >> (Array<T,2>& p){
  int ww = width_;
  for(int ii=0; ii < p.size1(); ii++){
    width_ = 2;
    put(ii);
    put(":");
    width_ = ww;
    for(int jj=0; jj < p.size2(); jj++)
      put(p(ii,jj));
    put("\n");
  }
  return *this;
}

template <class T, int Dim> Ostring& Ostring::operator >> (Tiny<T,Dim>& p){
  for(int ii=0; ii < Dim; ii++)
    put(p(ii));
  return *this;
}

#endif 
