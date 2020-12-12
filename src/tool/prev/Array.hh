#ifndef __ARRAY_MEC_HH_
#define __ARRAY_MEC_HH_

#define DIAG_ERROR 1

#include <cstdio>
#include <ctime>
#include <iostream>
#include <cmath>
#include <complex>
#include "Basic.hh"
//#include "Diag.hh"
#include "String.hh"
#include "Token.hh"
#include "Factory.hh"
#include "Serial.hh"


template <class T>
T sqrt(T);
template <class T, int Dim>
class Array;
template <class T >
class Array<T,2>;

template <class V>
V norm_2(Array<V,1>& xvec,Array<V,1>& yvec);

template <class V>
void unit(Array<V,1>& xvec,  V& norm, Array<V,1>& yvec);

template <class V>
void identity(Array<V,2>& matrix);

//
//*********************************************      scale V V
//
template <class V>
int scale(V alpha, Array<V,1>& xvec, V beta,
                 Array<V,1>& yvec);
template <class V>
int scale(V alpha, Array<V,1>& xvec,
                 Array<V,1>& yvec);
template <class V>
int scale(V alpha, Array<V,1>& xvec);

template <class V>
int scale(       Array<V,1>& xvec,
                 Array<V,1>& yvec);
//
//*********************************************      scale M M
//

template <class V>
int scale(V alpha, Array<V,2>& xvec, V beta,
                 Array<V,2>& yvec);
template <class V>
int scale(       Array<V,2>& xvec, V beta,
                 Array<V,2>& yvec);
template <class V>
int scale(V alpha, Array<V,2>& xvec);

template <class V>
int scale(V alpha, Array<V,2>& xvec,
                 Array<V,2>& yvec);
template <class V>
int scale(       Array<V,2>& xvec,
                 Array<V,2>& yvec);

//
//*********************************************      copy vec vec
//
template <class V, class T>
int copy(Array<V,1>& b_vec, int bib, int bie,
	 Array<T,1>& c_ve);
template <class V, class T>
int copy(Array<V,1>& b_vec, int bib, int bie,
	 Array<T,1>& c_ve, int cib);
template <class V, class T>
int copy(Array<V,1>& b_vec,
	 Array<T,1>& c_ve );
template <class V, class T>
int copy(Array<V,1>& b_vec,
	 Array<T,1>& c_ve, int cib);
template <class V, class T>
int copy(Array<V,1>& b_vec, int bib,
	 Array<T,1>& c_ve);
template <class V, class T>
int copy(Array<V,1>& b_vec, int bib,
	 Array<T,1>& c_ve, int cib);
//
//*********************************************      copy mat vec
//

template <class V, class T>
int copy(Array<V,2>& a_mat, int row,  Array<T,1>& b_vec);      //useful as copy from storage to vector, process

template <class V, class T>
int copy(Array<V,2>& a_mat, int row,int ajb, int aje,  Array<T,1>& b_vec);      //useful as copy from storage to vector, process

template <class V, class T>
int copy(Array<V,2>& a_mat, int aib, int aie, int ajb, int aje,
                     Array<T,1>& b_vec);


template <class V, class T>
int copy(Array<V,2>& a_mat, int aib, int aie, int ajb, int aje,
                     Array<T,1>& b_vec, int bib, int bie);

template <class V, class T>
int copy(Array<V,1>& b_vec,  Array<T,2>& a_mat, int row);      //useful as copy to storage from vector, process

template <class V, class T>
int copy(Array<V,1>& b_vec,  Array<T,2>& a_mat, int row, int ajb, int aje);      //useful as copy to storage from vector, process


template <class V, class T>
int copy(Array<V,1>& b_vec, 
                     Array<T,2>& a_mat, int aib, int aie, int ajb, int aje);

template <class V, class T>
int copy(Array<V,1>& b_vec, int bib, int bie,
                     Array<T,2>& a_mat, int aib, int aie, int ajb, int aje);

//
//*********************************************      copy mat mat
//
template <class V, class T>
int copy( Array<V,2>& xmat,Array<T,2>& ymat);

template <class V, class T>
int copy( Array<V,2>& xmat, 
          Array<T,2>& ymat, int iybeg);

template <class V, class T>
int copy( Array<V,2>& xmat, int ixbeg,
          Array<T,2>& ymat);


template <class V, class T>
int copy( Array<V,2>& xmat, int nrow, int ixbeg, int ixend,  //useful for storage
          Array<T,2>& ymat);

template <class V, class T>
int copy( Array<V,2>& xmat, int ixbeg, int ixend,
          Array<T,2>& ymat);

template <class V, class T>
int copy( Array<V,2>& xmat, int ixbeg, int iexend,
          Array<T,2>& ymat, int iybeg);

template <class V, class T>                                   // useful for storage
int copy( Array<V,2>& xmat, 
          Array<T,2>& ymat, int nrow, int iybeg, int iyend);

template <class V, class T>
int copy( Array<V,2>& xmat, 
          Array<T,2>& ymat, int iybeg, int iyend);

template <class V, class T>
int copy( Array<V,2>& xmat,
          Array<T,2>& ymat, int iybeg, int iyend, int jybeg, int jyend);

template <class V, class T>
int copy( Array<V,2>& xmat, int ixbeg, int ixend, int jxbeg, int jxend,
          Array<T,2>& ymat);

template <class V, class T>
int copy( Array<V,2>& xmat, int ixbeg, int ixend, int jxbeg, int jxend,
          Array<T,2>& ymat, int iybeg, int iyend, int jybeg, int jyend);

//
//*********************************************      
//

template <class V>
int prod(V alpha, Array<V,2>& mat, Array<V,1>& xvec,
                V beta,Array<V,1>& yvec);
template <class V>
int prod(V alpha, Array<V,1>& xvec, Array<V,2>& mat,
                V beta,Array<V,1>& yvec);
template <class V>
int prod(V alpha, Array<V,2>& mat, Array<V,2>& xvec,
                V beta,Array<V,2>& yvec);

template <class V>
int transpose(Array<V,2>& xmat, Array<V,2>& ymat);

template <class V>
int cross_prod(Array<V,1>& xvec,Array<V,1>& yvec,
                        Array<V,1>& zvec);

template <class V>
V dot_prod(Array<V,1>& xvec,Array<V,1>& yvec);

template <class V>
V dot_prod(Array<V,2>& xmat,Array<V,2>& ymat);

template <class V>
void set_zero(V& var);


//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//
//               -----  class Array  -----
//
//					    Copyright
//                                          Mechanics Computations,inc.
//					    Year 2000
//                                          All Rights Reserved
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

template <class T, int Dim>
class Array{

private:
  
  Array(Array<T,Dim>& ar){
    if(this != &ar){
      dim1 = ar.size();
      
      if(dim1 > maxrow){
	maxrow = dim1;
	if(buf ) delete[] buf;
	buf = new T [maxrow];
      }

      for(int i=0; i<maxrow; i++)
	buf[i] = ar.buf[i];
    }
  }

  Array& operator = (Array<T,Dim>& ar){

    if(this != &ar){
      dim1 = ar.size();

      if(dim1 > maxrow){
	maxrow = dim1;
	if(buf > 0 )   delete[] buf;
	buf = new T [maxrow];
      }

      for(int ii=0; ii<maxrow; ii++)
        buf[ii] = ar.buf[ii];
    }
    return *this;
  }

protected:

  friend class Analyzer;
  friend class Class;
  T *buf;
  int  dim1;
  int maxrow;

public:
  static bool is_initialize;
//
//************************************************      creation,copy,etc
//
  
  Array():dim1(0),maxrow(0),buf(0){
  }
  
  Array(int is1):dim1(is1),maxrow(is1),buf(0){
    if(maxrow > 0)
      buf = new T[maxrow]; 
    clear();
  }

  void dim(int dim_1)
  {
    dim1 = dim_1;

    if(dim1 > maxrow){
      maxrow = dim1;
      if(buf ) delete[] buf;
      buf = new T[maxrow]; 
    }
  
    clear();
  }
 
  ~Array(){
      if(buf) delete[] buf;
  }

  
  void serial(){        // commented out for debugging leak
    serializer<dim1,maxrow;
    serializer.pointer(buf,maxrow);
  }
  int get_serial_id(){return TypeFactory<Array>::get_id();}

  
  Array& operator = (T ar){
    for(int ii=0; ii<dim1; ii++)
      buf[ii] = ar;
    return *this;
  }

  T* pointer(){
    return buf;
  }

  //  This makes the above subroutine segfault. Don't know why. They key is the definition of one more 
  //  integer, index. Array cannot take it. Tiny seems to working in the Color class. It was done way after this.
  //  Look at index initialization.
//@  Array& operator , (T val){
//@    if(index < dim1)  
//@      buf[index++] = val;
//@    else
//@      diag_mesg(DIAG_ERROR,"Array<T,dim> , index out of bound:"<<index<<", Dim:"<<dim1<<std::endl);
//@
//@    return *this;
//@  }
//
//************************************************      
//
  
  T& operator()(int i){
    if(i>=0 && i<dim1) return buf[i];
    else{
      diag_mesg(DIAG_ERROR,"Array:() index out of bounds:"<<i<<":dim1:"<<dim1<<std::endl);
      return buf[0];
    }
  }

//
//*********************************************      
//
  T norm(){
    T norm_val(0);
    for(int i=0; i<dim1; i++)
      norm_val += buf[i]*buf[i];
    norm_val = sqrt(norm_val);
    return norm_val;
  }
  
  void clear(){
    for(int i=0; i < maxrow; i++)
      set_zero(buf[i]);
  }
  
  int size(){return dim1;}
  int maxsize(){return maxrow;}
//
//*********************************************      see also redim
//
  void resize(int dim1_){

    if(dim1_ > maxrow){

      int maxrow_o = maxrow;
      maxrow       = std::max(dim1_,2*maxrow);
      T* aux        = new T[maxrow];


      for(int ii=0; ii<maxrow; ii++)
	set_zero(aux[ii]);

      for(int ii =0; ii < maxrow_o; ii++)
        aux[ii] = buf[ii];
      
      if(buf) delete[] buf;
      buf  = aux;

    }
    
    dim1 = dim1_;
    
  }
  
//
//*********************************************     dim never gets smaller, but checks for enough memory. good for analyzer
//
  void redim(int dim1_){
    dim1_ = std::max(dim1,dim1_);
    resize(dim1_);
  }

//
//*********************************************      
//
  
  void scan(String st, int begin = 0, int end = i_bign){
    Token tk(st);
    int i = 0;
    end   = std::min(dim1-1,end);
    begin = std::max(0,begin);

    for(int ii = begin; ii <= end; ii++)
      set_zero(buf[ii]);
    
    while(tk.next() & tk.number){
      if( i >= begin && i <= end)
        buf[i++] = tk.fvalue;
    }
    
  }


  
  operator const T*() const{
    return (const T*) buf;
  }
  
//
//************************************************      friends
//
  friend int cross_prod<T>(Array<T,1>& xvec,Array<T,1>& yvec,
                        Array<T,1>& zvec);
  
  friend int cross_prod(Array<double,1>& xvec,Array<Complex,1>& yvec,
  			Array<Complex,1>& zvec);

  friend  T dot_prod<T>(Array<T,1>& xvec,Array<T,1>& yvec);

  friend int scale<T>(T alpha, Array<T,1>& xvec, T beta,
                   Array<T,1>& yvec);

  friend T norm_2<T>(Array<T,1>& xvec,Array<T,1>& yvec);
  
  friend void unit<T>(Array<T,1>& xvec, T& norm, Array<T,1>& yvec);
  
  friend int prod<T>(T alpha, Array<T,1>& xvec, Array<T,2>& mat, 
                  T beta,Array<T,1>& yvec);
  friend int prod<T>(T alpha, Array<T,2>& mat, Array<T,1>& xvec,
                  T beta,Array<T,1>& yvec);

  friend void set_zero<T>(T& var);
//
//************************************************      
//
  static bool initialize(){
    return Factory< Array<T,Dim> >::initialize();
  }
//
//************************************************      
//
};

template <class T, int Dim>
bool Array<T,Dim>::is_initialize = Array<T,Dim>::initialize();

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//
// -----  template <class T, int dim1>class Tiny:public Array<T,1>  -----
//
//					    Copyright
//                                          Mechanics Computations,inc.
//					    Year 2000
//                                          All Rights Reserved
//
//
// C: Tiny is a fix dim 1-d Array. Very handy for coords. Inherits all above.
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

template <class T, int Dim>
class Tiny:public Array<T,1>{

private:

  Tiny(const Tiny<T,Dim>& ar):Array<T,1>(Dim){
      Array<T,1>::index = Dim;
    for(int ii=0; ii<Dim; ii++)
      this->buf[ii] = ar.buf[ii];
  }


public:
  static bool is_initialize;

  Tiny():Array<T,1>(Dim){
  }
    
  Tiny& operator = (Tiny<T,Dim>& ar){

    if(this != &ar ){
      for(int ii=0; ii<Dim; ii++)
	this->buf[ii] = ar.buf[ii];
    }
    return *this;
  }

  Tiny& operator = (T val){
    for(int ii=0; ii<Dim; ii++)
      this->buf[ii] = val;
    return *this;
    }
  

  void serial(){ // commented out for debugging leak
    serializer.pointer(this->buf,this->dim1);
  }
  int get_serial_id(){return TypeFactory<Tiny>::get_id();}
  
//
//************************************************      
//
  static bool initialize(){
    return Factory< Tiny<T,Dim> >::initialize();
  }

};

template <class T, int dim1>
bool Tiny<T,dim1>::is_initialize = Tiny<T,dim1>::initialize();

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//
//               -----  template < class T >class Array<T,2>  -----
//
//					    Copyright
//                                          Mechanics Computations,inc.
//					    Year 2000
//                                          All Rights Reserved
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

template < class T >
class Array<T,2>{
  
protected:
  

  friend class Class;

  T *buf;
  int  dim1,dim2;
  int maxrow,maxcol;

public:

  Array(int dim1_, int dim2_):
    dim1(dim1_),dim2(dim2_),maxrow(dim1_),maxcol(dim2_),buf(0){

    if(maxrow*maxcol > 0)
      buf = new T [maxrow*maxcol];
    clear();
  }
  
  Array():
    dim1(0),dim2(0),maxrow(0),maxcol(0),buf(0){
    
  }

  ~Array(){
    if(buf) delete[] buf;
  }

  void dim(int dim_1, int dim_2)     // useful for static variables;
  {
    dim1 = dim_1;
    dim2 = dim_2;

    if(dim1 > maxrow || dim2 > maxcol){
      if(buf) delete[] buf;
      maxrow =  dim1;
      maxcol  = dim2;
      buf     = new T[maxrow*maxcol]; 
    }
    clear();
  }

  T* pointer(){
    return buf;
  }
  
  void serial(){ // commented out for debugging leak
    serializer<dim1,dim2,maxrow,maxcol;
    serializer.pointer(buf,maxrow*maxcol);
  }
  int get_serial_id(){return TypeFactory<Array>::get_id();}

  
  int size1(){return dim1;}
  int size2(){return dim2;}
  int size(){return dim1*dim2;}

  int maxsize1(){return maxrow;}
  int maxsize2(){return maxcol;}
  int maxsize(){return  maxrow*maxcol;}
  
//
//*********************************************      see also redim
//
  void resize(int dim1_){                      // useful for buffers that contain fix size vectors.
    resize(dim1_,dim2);
  }

  void resize(int dim1_,int dim2_){

    bool is_resize = false;
    int maxcol_o = maxcol;
    int maxrow_o = maxrow;
    
    T* aux;

    if(dim1_ > maxrow ){ 
      maxrow = std::max(dim1_,2*maxrow);
      is_resize = true;
    }
    if(dim2_ > maxcol ) {
      maxcol = std::max(dim2_,2*maxcol);
      is_resize = true;
    }
    
    if(is_resize){

      aux  = new T[maxrow*maxcol];                // always redo: may change the relative position of elements

      for(int ii=0; ii<maxrow*maxcol; ii++)
        set_zero(aux[ii]);
    
      for(int ii =0; ii < maxrow_o; ii++)
        for(int jj =0; jj < maxcol_o; jj++)
          aux[ii*maxcol+jj] = buf[ii*maxcol_o+jj];
      
      if(buf) delete[] buf;
      buf  = aux;
    }

    dim1 = dim1_;
    dim2 = dim2_;
    
  }

//
//*********************************************      redim: the size never gets smaller.
//

  T& operator()(int i, int j){
    if(dim2*i+j < dim1*dim2)
      return buf[maxcol*i+j];
    else{
      diag_mesg(DIAG_ERROR,"Array<T,2>:(i,j) Out of bounds:"<<i<<":"<<j<<std::endl);
      return buf[0];
    }
  }
  

  void clear(){
    for(int i=0; i < maxrow*maxcol; i++)
      set_zero(buf[i]);
  }
  
  Array<T,2>& operator = (T ar){
    for(int i=0; i<dim1; i++)
      for(int j=0; j<dim2; j++)
        buf[i*maxcol+j] = ar;
    return *this;
  }
  
  operator const T*() const{
    return (const T*) buf;
  }
//
//************************************************      friend
//
  friend int scale<T>(T alpha, Array<T,2>& xmat, T beta,
                   Array<T,2>& ymat);

  friend  T dot_prod<T>(Array<T,2>& xmat,Array<T,2>& ymat);

  friend int prod<T>(T alpha,
                  Array<T,2>& a_mat,
                  Array<T,2>& b_mat,
                  T beta,
                  Array<T,2>& c_mat);

  friend int transpose<T>(Array<T,2>& xmat, Array<T,2>& ymat);

  friend int prod<T>(T alpha, Array<T,2>& mat, Array<T,1>& xvec,
                  T beta,Array<T,1>& yvec);

  friend int prod<T>(T alpha, Array<T,1>& xvec, Array<T,2>& mat, 
                  T beta,Array<T,1>& yvec);

  friend void set_zero<T>(T& var);

//
//************************************************      serial
//
  static bool initialize(){
    return Factory< Array<T,2> >::initialize();
  }
};

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//
// -----  template <class T, int dim1>class Tiny:public Array<T,1>  -----
//
//					    Copyright
//                                          Mechanics Computations,inc.
//					    Year 2000
//                                          All Rights Reserved
//
//
// C: Tiny is a fix dim 1-d Array. Very handy for coords. Inherits all above.
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

template <class T, int Dim1, int Dim2>
class Tiny2D:public Array<T,2>{

private:


  Tiny2D(const Tiny2D<T,Dim1,Dim2>& ar):Array<T,2>(Dim1,Dim2){

    for(int ii=0; ii<Dim1*Dim2; ii++)
      this->buf[ii] = ar.buf[ii];
  }

public:

  Tiny2D():Array<T,2>(Dim1,Dim2){
  }

  Tiny2D& operator = (T val){
    for(int ii=0; ii<Dim1*Dim2; ii++)
      this->buf[ii] = val;
    return *this;
  }

  Tiny2D& operator = (Tiny2D<T,Dim1,Dim2>& ar){
    if(this != &ar ){
      for(int ii=0; ii<Dim1*Dim2; ii++)
	this->buf[ii] = ar.buf[ii];
    }
    return *this;
  }


};

typedef Tiny<double,2>     Coor2  , A2Double;
typedef Tiny<int,2>        ICoor2 , A2Int;
typedef Tiny2D<double,2,2> M2Double;
typedef Tiny2D<int,2,2>    M2Int;
typedef Tiny2D<Complex,2,2>M2Complex;
typedef Array< Coor2,1>    ACoor2 ;
typedef Array<ICoor2,1>    AICoor2;
typedef Array< Coor2,2>    MCoor2 ;
typedef Array<ICoor2,2>    MICoor2;

typedef Tiny<double,3>     Coord,  A3Double  ;
typedef Tiny<Complex,3>    CCoord, A3Complex  ;
typedef Tiny<int,3>        ICoord, A3Int ;
typedef Tiny<bool,3>       BCoord;
typedef Tiny2D<double,3,3> M3Double;
typedef Tiny2D<int,3,3>    M3Int;
typedef Tiny2D<Complex,3,3>M3Complex;
typedef Array< Coord,1>    ACoord ;
typedef Array<ICoord,1>    AICoord;
typedef Array< Coord,2>    MCoord ;
typedef Array<ICoord,2>    MICoord;

typedef Tiny<double,4>      Coor4,  A4Double;
typedef Tiny<int,4>         ICoor4, A4Int;
typedef Tiny<Complex,4>     CCoor4, A4Complex;
typedef Tiny2D<double,4,4>  M4Double;
typedef Tiny2D<int,4,4>     M4Int;
typedef Tiny2D<Complex,4,4> M4Complex;
typedef Array< Coor4,1>     ACoor4 ;
typedef Array<ICoor4,1>     AICoor4;
typedef Array< Coor4,2>     MCoor4 ;
typedef Array<ICoor4,2>     MICoor4;

typedef Tiny<double,6>      Coor6,  A6Double;
typedef Tiny<int,6>         ICoor6, A6Int;
typedef Tiny<bool,6>        BCoor6;
typedef Tiny<Complex,6>     CCoor6, A6Complex;
typedef Tiny2D<double,6,6>  M6Double;
typedef Tiny2D<int,6,6>     M6Int;
typedef Tiny2D<Complex,6,6> M6Complex;
typedef Array< Coor6,1>     ACoor6 ;
typedef Array<ICoor6,1>     AICoor6;
typedef Array< Coor6,2>     MCoor6 ;
typedef Array<ICoor6,2>     MICoor6;

typedef Tiny2D<double,12,12> M12Double;

typedef Array<bool,1>    ABool;
typedef Array<int,1>     AInt;
typedef Array<void*,1>   APointer;
typedef Array<double,1>  ADouble;
typedef Array<Complex,1> AComplex;
typedef Array<String,1>  AString;

typedef Array<bool,2>    MBool;
typedef Array<int,2>     MInt;
typedef Array<double,2>  MDouble;
typedef Array<Complex,2> MComplex;
typedef Array<String,2>  MString;


//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//
//
//					    Copyright
//                                          Mechanics Computations,inc.
//					    Year 2000
//                                          All Rights Reserved
//
//
// C: 
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

//
//
//************************************************      
//


template <class V>
int cross_prod(Array<V,1>& xvec,Array<V,1>& yvec,
          Array<V,1>& zvec){
//
//************************************************      compute, NO OVERWRITE.
//
  Array<V,1> aux(3);
  if(xvec.size() == 3 && yvec.size() == 3 && zvec.size() == 3){
    aux.buf[0] = xvec.buf[1]*yvec.buf[2] - xvec.buf[2]*yvec.buf[1]; 
    aux.buf[1] = xvec.buf[2]*yvec.buf[0] - xvec.buf[0]*yvec.buf[2]; 
    aux.buf[2] = xvec.buf[0]*yvec.buf[1] - xvec.buf[1]*yvec.buf[0];
    scale(aux,zvec);
  }
  else{
    diag_mesg(DIAG_ERROR,"Cross product of vectors with dim neq 3:"
              <<xvec.size()<<","<<yvec.size()<<","<<zvec.size()<<std::endl);
    return 0;
  }
//
//************************************************      done
//
  return 1;
}


//
//************************************************      dot product
//
template <class V>
V dot_prod(Array<V,1>& xvec,Array<V,1>& yvec){
//
//************************************************      set
//
  V dot(0);
//
//************************************************      compute
//
  if(xvec.dim1 == yvec.dim1){
    for(int i=0; i<xvec.dim1; i++)
      dot += xvec.buf[i]*yvec.buf[i];
  }
  else
    diag_mesg(DIAG_ERROR,"Dot product of vectors with different size\n");
//
//************************************************      done
//
  return dot;
}

template <class V>
V dot_prod(Array<V,2>& xmat,Array<V,2>& ymat){
//
//************************************************      set
//
  V dot(0);
//
//************************************************      compute
//
  if(xmat.dim1 == ymat.dim1 && xmat.dim2 == ymat.dim2 ){
    for(int ii=0; ii<xmat.dim1; ii++)
      for(int jj=0; jj<xmat.dim2; jj++)
	dot += xmat.buf[ii*xmat.maxcol+jj]*ymat.buf[ii*ymat.maxcol+jj];
  }
  else
    diag_mesg(DIAG_ERROR,"Dot product of matrices with different size\n");
//
//************************************************      done
//
  return dot;
}

template <> inline
Complex dot_prod(Array<Complex,2>& xmat,Array<Complex,2>& ymat){
//
//************************************************      set
//
  Complex dot(0);
//
//************************************************      compute
//
  if(xmat.dim1 == ymat.dim1 && xmat.dim2 == ymat.dim2 ){
    for(int ii=0; ii<xmat.dim1; ii++)
      for(int jj=0; jj<xmat.dim2; jj++)
	dot += xmat.buf[ii*xmat.maxcol+jj] * std::conj(ymat.buf[ii*ymat.maxcol+jj]);
  }
  else
    diag_mesg(DIAG_ERROR,"Dot product of matrices with different size\n");
//
//************************************************      done
//
  return dot;
}
//
//************************************************      dot product
//

template <class V>
V norm_2(Array<V,1>& xvec,Array<V,1>& yvec){
//
//************************************************      set
//
  V dot(0);
//
//************************************************      compute
//
  if(xvec.dim1 == yvec.dim1){
    for(int i=0; i<xvec.dim1; i++)
      dot += (xvec.buf[i] - yvec.buf[i])*(xvec.buf[i] - yvec.buf[i]);
  }
  else
    diag_mesg(DIAG_ERROR,"Distance of vectors with different size\n");
//
//************************************************      done
//
  return sqrt(dot);
}


template <class V>

void  unit(Array<V,1>& xvec,       V& vnorm, Array<V,1>& yvec){
//
//************************************************      set
//
  if(xvec.dim1 == yvec.dim1){
    set_zero(vnorm);
    vnorm = dot_prod(xvec,xvec);
    if(vnorm < d_eps){
      diag_mesg(diag.detail,"Unit vector: vector components almost zero\n");
      vnorm = 0.0;
      yvec  = 0.e0;
    }
    else{
      vnorm = sqrt(vnorm);
      for(int i=0; i<xvec.dim1; i++)
        yvec.buf[i] = xvec.buf[i]/vnorm;
    }
  }
  else
    diag_mesg(DIAG_ERROR,"Distance of vectors with different size\n");
//
//************************************************      done
//
}

template <class V>
void identity(Array<V,2>& matrix)
{
  if(matrix.size1() != matrix.size2()){
    diag_mesg(DIAG_ERROR,"Identity: matrix is not square:"<<matrix.size1()<<" : "<<matrix.size2()<<std::endl);
    return;
  }
  matrix = d_zero;
  for(int ii=0; ii<matrix.size1(); ii++)
    matrix(ii,ii) = d_one;

}
//
//************************************************      y = alpha x + beta Y
//
template <class V>
int scale(Array<V,1>& xvec,
          Array<V,1>& yvec){
  
  V alpha(1),beta(0);
  return scale(alpha,xvec,beta,yvec);
}

template <class V>
int scale(V alpha , Array<V,1>& xvec,
          Array<V,1>& yvec){
  
  V beta(0);
  return scale(alpha,xvec,beta,yvec);
  
}

template <class V>
int scale(V alpha , Array<V,1>& xvec)
{
  V beta(0);
  return scale(alpha,xvec,beta,xvec);
  
}

template <class V>
int scale(Array<V,1>& xvec, V beta,
          Array<V,1>& yvec){
  
  V alpha(1);
  return scale(alpha,xvec,beta,yvec);
  
}

template <class V>
int scale(V alpha, Array<V,1>& xvec, V beta,
          Array<V,1>& yvec){
  V aux;

  if(yvec.dim1 == 0)                      // interpret: on purpose
    yvec.resize(xvec.dim1);

  if( xvec.dim1 == yvec.dim1){
    for(int i = 0; i< xvec.dim1; i++){
      aux =  std::abs(beta) ?  beta*yvec.buf[i] : 0 ;
      yvec.buf[i] = aux + alpha*xvec.buf[i];
    }
    return 0;
  }
  else{
    diag_mesg(DIAG_ERROR,"Array: Trying to scale two vectors with different size\n");
    return 1;
  }
}

template <> inline
int scale( Array<String,1>& xvec,
          Array<String,1>& yvec){

  if( xvec.size() == yvec.size()){
    for(int ii=0; ii<xvec.size(); ii++)
      yvec(ii) = xvec(ii);
    return 0;
  }
  else{
    diag_mesg(DIAG_ERROR,"Array: Trying to scale two vectors with different size\n");
    return 1;
  }
}

//
//************************************************      scale <V,2>
//
template <class V>
int scale(Array<V,2>& xmat,
          Array<V,2>& ymat){
  
  V alpha(1),beta(0);
  return scale(alpha,xmat,beta,ymat);
}

template <class V>
int scale(V alpha , Array<V,2>& xmat,
          Array<V,2>& ymat){
  
  V beta(0);
  return scale(alpha,xmat,beta,ymat);
  
}

template <class V>
int scale(V alpha , Array<V,2>& xmat){
  V beta(0);
  return scale(alpha,xmat,beta,xmat);
}


template <class V>
int scale(Array<V,2>& xmat, V beta,
          Array<V,2>& ymat){
  
  V alpha(1);
  return scale(alpha,xmat,beta,ymat);
  
}

//
//*********************************************      x and y can be the same 
//

template <class V>
int scale(V alpha, Array<V,2>& xmat, V beta,
          Array<V,2>& ymat){

  int dim1 = xmat.dim1, dim2=xmat.dim2;
  V aux;

  if( dim1 == ymat.dim1 &&  dim2 == ymat.dim2 ){
    for(int ii = 0; ii< dim1; ii++)
      for(int jj = 0; jj< dim2; jj++){
	aux  =  std::abs(beta) ?  beta*ymat.buf[ii*ymat.maxcol+jj] : 0 ;
	ymat.buf[ii*ymat.maxcol+jj] = aux + alpha*xmat.buf[ii*xmat.maxcol+jj];
      }
    return 0;
  }
  else{
    diag_mesg(DIAG_ERROR,"Array: Trying to scale two matrices with different size\n");
    return 1;
  }
}

//
//*********************************************      copy v-v
//
template <class V, class T>
int copy(Array<V,1>& bvec,
	 Array<T,1>& cvec)
{
  return copy(bvec,0,bvec.size()-1,cvec,0);
}
template <class V, class T>
int copy(Array<V,1>& bvec, 
	 Array<T,1>& cvec, int cib)
{
  return copy(bvec,0,bvec.size()-1,cvec,cib);
}
template <class V, class T>
int copy(Array<V,1>& bvec, int bib,
	 Array<T,1>& cvec)
{
  return copy(bvec,bib,bvec.size()-1,cvec,0);
}
template <class V, class T>
int copy(Array<V,1>& bvec, int bib, int bie,
	 Array<T,1>& cvec)
{
  return copy(bvec,bib,bie,cvec,0);
}
template <class V,class T>
int copy(Array<V,1>& bvec, int bib,
	 Array<T,1>& cvec, int cib)
{
  return copy(bvec,bib,bvec.size()-1,cvec,cib);
}
template <class V, class T>
int copy(Array<V,1>& bvec, int bib, int bie,
	 Array<T,1>& cvec, int cib)
{
  int num = bie-bib+1;
  int blast = bib+num;
  int clast = cib+num;
  if(blast > bvec.size() || clast > cvec.size()){
    diag_mesg(DIAG_ERROR,"Array::copy: Index out of bounds: bvec:"<<bvec.size()<<"    "<<cvec.size()<<std::endl);
    return 0;
  }
  for(int ii=0; ii< bie-bib+1; ii++)
    cvec(ii+cib) = bvec(ii+bib);
  return 1;
}
//
//*********************************************      copy mat vec
//
template <class V, class T>
int copy(Array<V,2>& a_mat, int row,  Array<T,1>& b_vec)      //useful as copy from storage to vector, process
{
  int acol = b_vec.size()-1;
  int bcol = a_mat.size2()-1;
  return copy(a_mat,row,row,0,acol, b_vec,0,bcol);
}

template <class V, class T>
int copy(Array<V,2>& a_mat, int row,  int ajb, int aje, Array<T,1>& b_vec)          //useful as copy from storage to vector, process
{
  return copy(a_mat,row,row,ajb,aje,  b_vec,0,b_vec.size()-1);
}


template <class V, class T>
int copy(Array<V,2>& a_mat, int aib, int aie, int ajb, int aje,
                     Array<T,1>& b_vec)                            
{
  return copy(a_mat,aib,aie,ajb,aje,  b_vec,0,b_vec.size()-1);
}


template <class V, class T>
int copy( Array<V,2>& xmat, int ixbeg, int ixend,
                            int jxbeg, int jxend,
          Array<T,1>& yvec, int iybeg, int iyend){
  
  if(std::min(ixbeg,ixend)<0 || std::max(ixbeg,ixend) >= xmat.size1() ||
     std::min(jxbeg,jxend)<0 || std::max(jxbeg,jxend) >= xmat.size2() ||
     std::min(iybeg,iyend)<0 || std::max(iybeg,iyend) >= yvec.size() ){
    
    diag_mesg(DIAG_ERROR,"Array:copy:Index out of bounds: xmat:size1= "<<xmat.size1()<<" :"<<ixbeg<<","<<ixend<<std::endl<<
	      "                                xmat.size2= "<<xmat.size2()<<" :"<<jxbeg<<","<<jxend<<std::endl<<
	      "                                yvec.size= "<<yvec.size()<<" :"<<iybeg<<","<<iyend<<std::endl);
    return 0;
  }
  
  else if( (ixend-ixbeg+1)*(jxend-jxbeg+1) != (iyend-iybeg+1)){
    diag_mesg(DIAG_ERROR,"Array:copy:Number of objects differ:"<<
              (ixend-ixbeg+1)*(jxend-jxbeg+1)<<":"<<(iyend-iybeg+1));
    return 0;
  }
  
  else{
    int kk = 0;
    for(int ii = ixbeg; ii<= ixend; ii++)
      for(int jj = jxbeg; jj<= jxend; jj++)
        yvec(kk++) = xmat(ii,jj) ;
  }
  
  return 1;
}

//
//*********************************************      copy vec mat
//
template <class V, class T>
int copy(Array<V,1>& b_vec, Array<T,2>& a_mat, int row)      //useful as copy from storage to vector, process
{
  int acol = b_vec.size()-1;
  int bcol = a_mat.size2()-1;
  return copy(b_vec,0,bcol, a_mat,row,row,0,acol);
}

template <class V, class T>
int copy(Array<V,1>& b_vec, Array<T,2>& a_mat, int row, int ajb, int aje) //useful as copy from storage to vector, process
{
  return copy(b_vec,0,b_vec.size()-1,  a_mat,row,row,ajb,aje);
}

template <class V, class T>
int copy(Array<V,1>& b_vec, Array<T,2>& a_mat, int aib, int aie, int ajb, int aje)
{
  return copy(b_vec,0,b_vec.size()-1,  a_mat,aib,aie,ajb,aje);
}


template <class V, class T>
int copy( Array<V,1>& yvec, int iybeg, int iyend,
          Array<T,2>& xmat, int ixbeg, int ixend,
                            int jxbeg, int jxend){

  
  if(std::min(ixbeg,ixend)<0 || std::max(ixbeg,ixend) >= xmat.size1() ||
     std::min(jxbeg,jxend)<0 || std::max(jxbeg,jxend) >= xmat.size2() ||
     std::min(iybeg,iyend)<0 || std::max(iybeg,iyend) >= yvec.size() ){
    diag_mesg(DIAG_ERROR,"Array:copy:Index out of bounds:yvec.size= "<<yvec.size()<<" :"<<iybeg<<","<<iyend<<std::endl<<
	                 "	                       xmat:size1= "<<xmat.size1()<<" :"<<ixbeg<<","<<ixend<<std::endl<<
	      "                               xmat.size2= "<<xmat.size2()<<" :"<<jxbeg<<","<<jxend<<std::endl);
    return 0;
  }
  
  else if( (ixend-ixbeg+1)*(jxend-jxbeg+1) != (iyend-iybeg+1)){
    diag_mesg(DIAG_ERROR,"Array:copy:Number of objects differ:"<<
              (ixend-ixbeg+1)*(jxend-jxbeg+1)<<":"<<(iyend-iybeg+1));
    return 0;
  }
  
  else{
    int kk = 0;
    for(int ii = ixbeg; ii<= ixend; ii++)
      for(int jj = jxbeg; jj<= jxend; jj++){
         xmat(ii,jj) = yvec(kk++)  ;
      }
  }
  
  return 1;
}

//
//*********************************************      copy mat - mat
//


template <class V, class T>
int copy( Array<V,2>& xmat,Array<T,2>& ymat){
  return copy(xmat, 0,xmat.size1()-1,     0,xmat.size2()-1, ymat, 0,xmat.size1()-1,          0,xmat.size2()-1);
}

template <class V, class T>
int copy( Array<V,2>& xmat, 
          Array<T,2>& ymat, int iybeg){
  return copy(xmat, 0,xmat.size1()-1,     0,xmat.size2()-1, ymat, iybeg,iybeg+xmat.size1()-1, 0,xmat.size2()-1);
}

template <class V, class T>
int copy( Array<V,2>& xmat, int ixbeg,
          Array<T,2>& ymat){
  return copy(xmat, ixbeg,xmat.size1()-1, 0,xmat.size2()-1, ymat, 0,xmat.size1()-1-ixbeg,     0,xmat.size2()-1);
}


template <class V, class T>
int copy( Array<V,2>& xmat, int row, int iybeg, int iyend,
          Array<T,2>& ymat){
  return copy(xmat, row,row, iybeg,iyend,    ymat, 0, ymat.size1()-1, 0, ymat.size2()-1);
}

template <class V, class T>
int copy( Array<V,2>& xmat, int ixbeg, int ixend,
          Array<T,2>& ymat){
  return copy(xmat, ixbeg,ixend,          0,xmat.size2()-1, ymat, 0,ixend-ixbeg,             0,xmat.size2()-1);
}

template <class V, class T>
int copy( Array<V,2>& xmat, 
          Array<T,2>& ymat, int row, int iybeg, int iyend){
  return copy(xmat, 0,xmat.size1()-1, 0,xmat.size2()-1,    ymat, row, row, iybeg, iyend);
}

template <class V, class T>
int copy( Array<V,2>& xmat, 
          Array<T,2>& ymat, int iybeg, int iyend){
  return copy(xmat, 0,xmat.size1()-1, 0,xmat.size2()-1,    ymat, iybeg, iyend, 0,xmat.size2()-1);
}

template <class V, class T>
int copy( Array<V,2>& xmat, int ixbeg, int ixend,
          Array<T,2>& ymat, int iybeg){
  return copy(xmat, ixbeg,ixend,          0,xmat.size2()-1, ymat, iybeg,iybeg+ixend-ixbeg,             0,xmat.size2()-1);
}



template <class V, class T>
int copy( Array<V,2>& xmat,
          Array<T,2>& ymat, int iybeg, int iyend, int jybeg, int jyend){
  return copy(xmat, 0,iyend-iybeg, 0,jyend-jybeg,    ymat, iybeg, iyend, jybeg, jyend);
}

template <class V, class T>
int copy( Array<V,2>& xmat, int ixbeg, int ixend, int jxbeg, int jxend,
          Array<T,2>& ymat){
  return copy(xmat, ixbeg,ixend,          jxbeg,jxend, ymat, 0,ixend-ixbeg,             0,jxend-jxbeg);
}


template <class V, class T>
int copy( Array<V,2>& xmat, int ixbeg, int ixend, int jxbeg, int jxend,
          Array<T,2>& ymat, int iybeg, int iyend, int jybeg, int jyend){

  if(std::min(ixbeg,ixend)<0 || std::max(ixbeg,ixend) >= xmat.size1() ||
     std::min(jxbeg,jxend)<0 || std::max(jxbeg,jxend) >= xmat.size2() ||
     std::min(iybeg,iyend)<0 || std::max(iybeg,iyend) >= ymat.size1() ||
     std::min(jybeg,jyend)<0 || std::max(jybeg,jyend) >= ymat.size2()
    ){
    diag_mesg(DIAG_ERROR,"Array:copy:Index out of bounds:xmat:size1= "<<xmat.size1()<<" :"<<ixbeg<<","<<ixend<<std::endl<<
                         "                               xmat.size2= "<<xmat.size2()<<" :"<<jxbeg<<","<<jxend<<std::endl<<
                         "                               ymat:size1= "<<ymat.size1()<<" :"<<iybeg<<","<<iyend<<std::endl<<
                         "                               ymat.size2= "<<ymat.size2()<<" :"<<jybeg<<","<<jyend<<std::endl

    );
    return 0;
  }

  else if( (ixend-ixbeg+1)*(jxend-jxbeg+1) != (iyend-iybeg+1)*(jyend-jybeg+1)){

    diag_mesg(DIAG_ERROR,"Array:copy:Number of objects differ:"<<((ixend-ixbeg+1)*(jxend-jxbeg+1))
        <<"!="<<((iyend-iybeg+1)*(jyend-jybeg+1))<<std::endl);
    return 0;
  }

  else{
    int ncol = jyend-jybeg+1;

    int kk = 0;
    for(int ii = ixbeg; ii<= ixend; ii++)
      for(int jj = jxbeg; jj<= jxend; jj++){
	int irow = kk / ncol + iybeg;
	int icol = kk % ncol + jybeg;
        ymat(irow,icol) =  xmat(ii,jj)  ;
	kk++;
      }
  }

  return 1;
}
//
//************************************************      m-v product
//

template <class V>
int prod(Array<V,2>& A,
         Array<V,1>& x,
         Array<V,1>& y){
  
  V alpha(1),beta(0);
  prod(alpha,A,x,beta,y);
  
  return 1;
}

template <class V>
int prod(Array<V,2>& A,
         Array<V,1>& x,
         V beta,
         Array<V,1>& y){
  
  V alpha(1);
  prod(alpha,A,x,beta,y);
  
  return 1;
}


template <class V> int prod(V alpha,
         Array<V,2>& mat,
         Array<V,1>& xvec,
         V beta,
         Array<V,1>& yvec){
  
  if(mat.dim2 != xvec.dim1 || mat.dim1 != yvec.dim1){
    diag_mesg(DIAG_ERROR,
              "2D Array index and 1D array index are inconsistent"<<
              mat.dim2<<","<< xvec.dim1<<std::endl);
    return 1;
  }

  V aux;
  Array<V,1> aux_vec(yvec.dim1);

  for(int i = 0; i< mat.dim1; i++){
    aux = std::abs(beta) ?  beta*yvec.buf[i] : 0 ;
    for(int j = 0; j < mat.dim2; j++)
      aux += alpha*mat.buf[i*mat.maxcol+j]*xvec.buf[j];
    aux_vec.buf[i] = aux;
  }

  scale(aux_vec,yvec);

  return 0;
}
//
//************************************************      v-m product
//

template <class V>
int prod(Array<V,1>& x,
	 Array<V,2>& A,
         Array<V,1>& y){
  
  V alpha(1),beta(0);
  prod(alpha,x,A,beta,y);
  
  return 1;
}

template <class V>
int prod(Array<V,1>& x,
         Array<V,2>& A,
         V beta,
         Array<V,1>& y){
  
  V alpha(1);
  prod(alpha,x,A,beta,y);
  
  return 1;
}


template <class V> int prod(V alpha,
         Array<V,1>& xvec,
         Array<V,2>& mat,
         V beta,
         Array<V,1>& yvec){
  
  if(mat.dim1 != xvec.dim1 || mat.dim2 != yvec.dim1){
    diag_mesg(DIAG_ERROR,
              "2D Array index and 1D array index are inconsistent:"<<
              mat.dim1<<","<<xvec.dim1<<std::endl);
    return 1;
  }

  V aux;
  Array<V,1> aux_vec(yvec.dim1);


  for(int j = 0; j < mat.dim2; j++){
    aux = std::abs(beta) ?  beta*yvec.buf[j] : 0 ;
    for(int i = 0; i< mat.dim1; i++)
      aux += alpha*mat.buf[i*mat.maxcol+j]*xvec.buf[i];
    aux_vec.buf[j] = aux;
  }

  scale(aux_vec,yvec);

  return 0;
}
//
//************************************************      m-m product
//

template <class V>
int prod(Array<V,2>& A,
         Array<V,2>& B,
         Array<V,2>& C){
  
  V alpha(1),beta(0);
  prod(alpha,A,B,beta,C);
  
  return 1;
}

//
//*********************************************      a_mat can be b_mat or c_mat
//

template <class V>
int prod(V alpha,
         Array<V,2>& b_mat,
         Array<V,2>& c_mat,
         V beta,
         Array<V,2>& a_mat){
//
//************************************************      check indexes
//
  
  if( a_mat.dim1 == b_mat.dim1 && b_mat.dim2 == c_mat.dim1 &&
      a_mat.dim2 == c_mat.dim2){

    V aux;
    Array<V,2> aux_mat(a_mat.dim1,a_mat.dim2);

    for(int i = 0; i< a_mat.dim1; i++)
      for(int j = 0; j< a_mat.dim2; j++){
	aux = std::abs(beta) ?  beta*a_mat.buf[i*a_mat.maxcol+j] : 0 ;
	for(int k = 0; k < b_mat.dim2; k++)
	  aux += alpha * b_mat.buf[ i*b_mat.maxcol+k ] * c_mat.buf[ k*c_mat.maxcol+j] ;
	aux_mat.buf[i*a_mat.maxcol+j] = aux;
      }
    scale(aux_mat,a_mat);
  }
  else{
    diag_mesg(DIAG_ERROR,
    "2D Matrix-Matrix multiplication indexes are inconsistent\n");
    return 1;
  }

  return 0;
}

template <class V>
int transpose(Array<V,2>& a_mat,
	      Array<V,2>& b_mat)
{
//
//************************************************      check indexes
//
  
  if( a_mat.dim1 == b_mat.dim2 &&
      a_mat.dim2 == b_mat.dim1){

    Array<V,2> b_aux(b_mat.dim1,b_mat.dim2);
    for(int i = 0; i< a_mat.dim1; i++)
      for(int j = 0; j< a_mat.dim2; j++)
	b_aux.buf[j*b_mat.maxcol+i] = a_mat.buf[i*a_mat.maxcol+j];
    scale(b_aux,b_mat);
  }
  else{
    diag_mesg(DIAG_ERROR,
    "2D Matrix transpose indexes are inconsistent\n");
    return 1;
  }

  return 0;
}



//
//*********************************************      produces segs, see uoc::out instead
//
template <class V>
std::ostream& operator << (std::ostream& out, Array<V,2>& b_mat){
  for(int i=0; i<b_mat.size1(); i++){
    out<<i<<": ";
    for(int j=0; j<b_mat.size2(); j++){
      out<<b_mat(i,j);
      out<<", ";
    }
    out<<std::endl;
  }
  return out;
}

//
//*********************************************      produces segs, see uoc::out instead
//
template <class V, int Dim>
std::ostream& operator << (std::ostream& out, Tiny<V,Dim>& b_vec){
  int ii;
  out<<"<";

  for(ii=0; ii<b_vec.size(); ii++){
    out<<b_vec(ii);
    if(ii<b_vec.size()-1) out <<",";
  }
  out<<">";
  return out;
}
//
//*********************************************      produces segs, see uoc::out instead
//
template <class V>
std::ostream& operator << (std::ostream& out, Array<V,1>& b_vec){
  for(int i=0; i<b_vec.size(); i++){
    out<<b_vec(i);
    out<<", ";
  }
  return out;
}
//
//*********************************************      zero functions
//

template <class V>
void set_zero(V& var){
  var = 0;
}

template <> inline
void set_zero(double& var){
  var = 0.e0;
}

template <> inline
void set_zero(Complex& var){
  var = Complex(d_zero,d_zero);
}

template<> inline
void set_zero(String& var){
  var = s_zero;
}

template<> inline
void set_zero(bool& var){
  var = false;
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//               -----  template <class V> void new_array  -----
//
//
// C: A handy way to dimension c++ 2dimensional array from a pointer
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

template <class V> void new_array(V**& matrix, int nrow, int ncol)
{
  matrix    = new V*[nrow];      
  matrix[0] = new V [nrow*ncol];
  for(int ii=1; ii<nrow; ii++)
    matrix[ii] = matrix[ii-1] + ncol;
}

#endif






