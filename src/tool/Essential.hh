#ifndef __ESSENTIAL_MEC_HH__
#define __ESSENTIAL_MEC_HH__
//#pragma interface
//
//************************************************      includes
//
#include <string>
#include <climits>
#include <cfloat>
#include <iomanip> 
#include <complex>
//
//************************************************      includes
//
typedef int   db_int;
typedef long  db_long;
typedef std::complex<double> Complex;
//
//************************************************      buffer
//
#define  namdim  8
#define  nbuffer  1020
extern char str_buffer[nbuffer];
//
//**************************************************    type definitions
//
#ifdef cchar
#undef cchar
#endif
#define cchar (const char*)
#define Minmax(x,l,u) max(x,l) < (u) ? (x) : (u)
//
//************************************************      number char strings
//
#define numchs(opt) ( sizeof(opt)/sizeof(char*) )
//
//*********************************************      base file
//
//define __FULL_FUNCTION__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 __FUNCTION__ : __FILE__ __FUNCTION__ ) 
#define __FULL_FUNCTION__    strrchr(__FILE__, '/') + 1 << "::" << __FUNCTION__
//
//*********************************************      
//
#define u_printp(x)     std::clog<<" +:" <<__PRETTY_FUNCTION__<<" : "<< #x <<" : "<<std::clog.precision(6)<<std::fixed<<x<<std::endl
#define u_printp_h(x)   std::clog<< "***"<<__PRETTY_FUNCTION__<<" : "<< #x <<" : "<<std::endl<<std::clog.precision(6)<<std::fixed<<x<<std::endl
#define u_printp_t(x)   std::clog<< "***"<<__PRETTY_FUNCTION__<<" : "<< #x <<" : "<<std::endl<<std::endl
#define u_print(x)      std::clog<<" +:" <<__FULL_FUNCTION__<<"\t: "<< #x <<" : "<<std::clog.precision(4)<<std::fixed<<x<<std::endl
#define u_print_h(x)    std::clog<<" +:"<<__FULL_FUNCTION__<<"\t: "<< #x <<" : "<<std::endl<<std::clog.precision(4)<<std::fixed<<x<<std::endl
#define u_print_t(x)    std::clog<<" +:"<<__FULL_FUNCTION__<<"\t: "<< #x <<" : "<<std::endl<<std::endl
#define u_printf(x,d)   std::clog<<" +:" <<__FULL_FUNCTION__<<"\t: "<< #x <<" : "<<std::clog.precision(d)<<std::fixed<<x<<std::endl
#define u_printe(x,d)   std::clog<<" +:" <<__FULL_FUNCTION__<<"\t: "<< #x <<" : "<<std::clog.precision(d)<<std::scientific<<x<<std::endl
#define u_printf_h(x,d) std::clog<<" +:"<<__FULL_FUNCTION__<<"\t: "<< #x <<" : "<<std::endl<<std::clog.precision(d)<<std::fixed<<x<<std::endl
#define u_printe_h(x,d) std::clog<<" +:"<<__FULL_FUNCTION__<<"\t: "<< #x <<" : "<<std::endl<<std::clog.precision(d)<<std::scientific<<x<<std::endl
#define u_print_f       std::clog<<" +:"<<__PRETTY_FUNCTION__<<std::endl
#define u_sep           "    : "
#define u_diag(x)       std::clog<<u_sep<< #x <<" \t: "<<x<<std::endl
#define u_diag_h(x)     std::clog<<u_sep<< #x <<" \t: "<<std::endl<<x<<std::endl
#define u_diag_t(x)     std::clog<<u_sep<< #x <<" \t: "<<std::endl<<std::endl
#define u_diage(x,d)    std::clog<<u_sep<< #x <<" \t: "<<diag.precision(d)<<std::scientific<<x<<std::endl
#define u_diagf(x,d)    std::clog<<u_sep<< #x <<" \t: "<<std::clog.precision(d)<<std::fixed<<x<<std::endl
#define u_diagf_h(x,d)  std::clog<<u_sep<< #x <<" \t: "<<std::endl<<std::clog.precision(d)<<std::scientific<<x<<std::endl
#define u_diag_f        std::clog<<" +++: "<<__FULL_FUNCTION__<<" : "<<std::endl
//
//************************************************      
//
#define diag_mesg(x,y) { std::clog << y ; }     //else : if doesnt mix well with other if above

#define u_program_error diag_list(diag.program,diag<<"Program error::"<< __PRETTY_FUNCTION__ <<std::endl);
//
//*********************************************      
//

#define diag_head(x) ("//")

#define method_delim   __PRETTY_FUNCTION__ 
#define file_delim "*****" __FILE__
#define RXint    String::rxint
#define RXdouble String::rxdouble
#define RXwhite  String::rxwhite
#define RXsep    String::rxseparator
//
//************************************************      number consts
//
extern int    i_zero;
extern int    i_one ;
extern int    i_bign;
extern int    i_litn;
extern int    i_undef;
extern int    i_argc;
extern int    i_all;
extern float  r_zero;
extern float  r_one ;
extern float  r_bign;    // 3.40282e+38
extern float  r_litn;    // 1.17549e-38
extern float  r_eps ;    // 1.19209e-07
extern float  r_epslen ;    // 1.19209e-07
extern float  r_epsang ;    // 1.19209e-07
extern float  r_small;   // 1.e-4
extern double d_zero;
extern double d_one ;
extern double d_half;
extern double d_big_unit; // 7777, values around unit going ballistic
extern double d_bign;    // 1.79769e+308
extern double d_litn;    // 2.22507e-308
extern double d_eps ;    // 2.22045e-16
extern double d_small ;    // 1.e-8
extern double d_epslen;  //
extern double d_epsang;  //
extern double d_epsf;   // small force, (for reporting).
extern double d_epsll;  // small length, 1cm
extern double d_epsgl;  // min length to show with gl, 0.1mm
extern double d_undef;
extern std::complex<double> z_zero;
extern std::complex<double> z_one;

const char c_zero            = '\0';
const char c_ends            = '\0';
const char c_endl            = '\n';
const char c_space           = ' ' ;
const char c_verbatim        = '\'';
const char c_tab             = '\t';
const char c_param_sub       = '%' ;
const char c_slash_dir       = '/' ;        // "\" for dos
const char c_forward_slash   = '/' ;           //  "/";
const char c_command         = '&' ;
const char c_comment         = '$' ;
const char c_class           = '~' ;
const char c_escape          = '\\';
const char c_backward_slash  = '\\';           //  "/";
const char c_comma           = ',' ;
const char c_quote           = '\"';
const char c_apostrophe      = '\'';
const char c_underscore      = '_';
const char c_dot             = '.';
const char c_semicolon       = ';';
const char c_colon           = ':';
const char c_column          = ';';
const char c_greater         = '>';
const char c_less            = '<';
const char c_internal_entity = '!';
//
//************************************************      strings
//
extern const char*  s_backward_slash;           //  "\";
extern const char*  s_category;
extern const char*  s_class;
extern const char*  s_comment;
extern const char*  s_command;
extern const char*  s_default;
extern const char*  s_dir_here;            //   "./"
extern const char*  s_dot;                     // "."
extern const char*  s_endl;
extern const char*  s_forward_slash;           //  "/";
extern const char*  s_quote;               // "\""
extern const char*  s_selall;
extern const char*  s_slash_dir;           //  "/";
extern const char*  s_space;               // " "
extern const char*  s_symmetry;
extern const char*  s_underscore;
extern const char*  s_undef;
extern const char*  s_verbatim;
extern const char*  s_vertical_slash;           //  "|";
extern const char*  s_white    ;           //  " \t\n\r";
extern const char*  s_zero     ;           //  ""     ;
extern const char*  s_yes;
extern const char*  s_no;
extern const char*  s_none;                // "none";
extern const char*  s_true;
extern const char*  s_false;
extern const char*  s_semicolon;           
extern const char*  s_colon;
extern const char*  s_global;              // "::"
extern const char*  s_column;              // ;
extern const char*  s_greater;             // ">";
extern const char*  s_less;                // "<";
//
//************************************************       constants
//
const double d_pi      = 3.14159265358979323846;
//
//************************************************      done
//
extern char**   c_argv;
extern int      i_argc;
const bool     b_true = true;
//
//*********************************************      
//
#endif

