//  calls with complex cannot be used this way. use this,
// #define zposv zposv_         // solve complex-double symmetric post definite
// extern "C"{
//  #include "Blas.h"
//  #include "Lapack.h"
//  int zposv_(char*,int*,int*,std::complex<double>*,int*,std::complex<double>*, int*, int*);
// }
//
#define dgeqrf dgeqrf_       // general qr factorization
#define dgetrf dgetrf_       // general factor
#define dgetri dgetri_       // general inverse after factor
#define dgetrs dgetrs_       // general solver
#define dormqr dormqr_       // action of dgecrf on a vector
#define dorgqr dorgqr_       // explicit form for Q after dgecrf
#define dpotrf dpotrf_       // factor pos definite symm
#define dpotri dpotri_       // inverse pos definite symm after dpotrf
#define dpotrs dpotrs_       // solve  pos def symm after dpotrf
#define dsytrd dsytrd_       // fact eigenvalues/vectors symm 
#define dorgtr dorgtr_       // eigenvalues after dsytrd
#define dsytrf dsytrf_       // factorize symm indefinite
#define dsytrs dsytrs_       // solve symm indefinite
#define dsytri dsytri_       // inverse symm indefinite after dsytrf
#define dtrtri dtrtri_       // computes the inverse of a triangular matrix
#define dsteqr dsteqr_       // obtain eigenvectors/values
#define dsyev dsyev_         // obtain eigenvectors/values with one call.

//char* uplo="U",*sidel="L",*transn="N",*transt="T";

int dgeqrf(int*,int*,double[],int*,double[], double[],int*,int*);
int dgetrf(int*, int*, double[], int*, int[], int*);
int dgetri(int*, double[], int*, int[], double[], int*, int*);
int dgetrs(char*, int*, int*, double[], int*, int*, double[], int*, int*);
int dormqr(char*,char*,int*,int*,int*,double[],int*,double[],double[],
	   int*,double[],int*,int*);
int dorgqr(int*,int*,int*,double[],int*,double[], double[],int*,int*);

int dsytrd(char*,int*,double[],int*,double[],double[],double[],double[],
	   int*,int*);
int dorgtr(char*,int*,double[],int*,double[],double[],int*,int*);
int dpotrf(char*,int*,double[],int*,int*);
int dpotri(char*,int*,double[],int*,int*);
int dpotrs(char*,int*,int*,double[],int*,double[],int*,int*);
int dsteqr(char*,int*,double[],double[],double[],int*,double[],int*);
int dsyev(char*,char*,int*,double[],int*,double[],double[],int*,int*);

int dsytrf(char*,int*,double[],int*,long int[],double[],int*,int*);
int dsytrs(char*,int*,int*,double[],int*,long int[],double[],int*,int*);
int dsytri(char*,int*,double[],int*,long int[],double[],int*);
int dtrtri(char*,char*,int*,double[],int*,int*);



