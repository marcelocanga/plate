#define ddot ddot_
#define dcopy dcopy_
#define dscal dscal_
#define daxpy daxpy_
#define dsymv dsymv_
#define dgemv dgemv_
#define dgemm dgemm_

double ddot(int*,double[],int*,double[],int*);
int dcopy(int*,double[],int*,double[],int*);
int dscal(int*,double*,double[],int*);
int daxpy(int*,double*,double[],int*, double[],int*);
int dsymv(char*,int*,double*,double[],int*,
	     double[],int*,double*,double[],int*);
int dgemv(char*,int*,int*,double*,double[],int*,
	     double[],int*,double*,double[],int*);
int dgemm(char*,char*,int*,int*,int*,double*,double[],int*,double[],int*,
	  double*,double[],int*);
