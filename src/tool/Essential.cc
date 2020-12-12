//#pragma implementation
#include "Essential.hh"

//
//************************************************      global buffer
//
int    i_zero  = 0 ;
int    i_one   = 1 ;
int    i_bign  =  INT_MAX;
int    i_litn  =  INT_MIN;
int    i_undef =  -77;
int    i_all   = 0x7FFFFFFF;
float  r_zero  = 0.    ;      
float  r_one   = 1.    ;
float  r_bign  = FLT_MAX;           // 3.40282e+38
float  r_litn  = FLT_MIN;           // 1.17549e-38
float  r_eps   = FLT_EPSILON;       // 1.19209e-07
float  r_epslen   = FLT_EPSILON*1.e3;       // 1.19209e-04
float  r_epsang   = FLT_EPSILON*1.e3;       // 1.19209e-04
float  r_small = 1.e-4;
double d_zero  = 0.e0  ;
double d_one   = 1.e0  ;
double d_half  = 0.5e0 ;
double d_big_unit = 7777.7777; // 7777, values around unit going ballistic
double d_bign  = 0.01  * DBL_MAX;           // 1.79769e+308
double d_litn  = 100.0 * DBL_MIN;           // 2.22507e-308
double d_eps   = DBL_EPSILON;       // 2.22045e-16
double d_epslen   = DBL_EPSILON*1.e6;       // 2.22045e-10
double d_epsang   = DBL_EPSILON*1.e6;       // 2.22045e-10
double d_small  = 1.e-8;
double d_epsgl  = 1.e-4;          // 1mm
double d_epsf  =  1.e-4;            // small force
double d_epsll =  1.e-2;            // small length, 1cm
double d_undef = -77.77;

std::complex<double> z_zero(0.0,0.0);
std::complex<double> z_one (1.0,0.0);

char str_buffer[nbuffer];


int i_argc;
char** c_argv;

const char* s_undef         = "####";  // dont change for @ $
const char* null_str        = "";  //deprec
const char* s_slash_dir     = "/";
const char* s_underscore    = "_";
const char* s_command       = "&";
const char* s_forward_slash = "/";
const char* s_vertical_slash= "|";
const char* s_backward_slash= "\\";
const char* s_white         = " \t\n\r";
const char* s_zero          = "";
const char* s_comment       = "$";
const char* s_space         = " ";
const char* s_default       = "default";
const char* s_class         = "~";
const char* s_category      = ":";
const char* s_selall        = ".*";
const char* s_endl          = "\n";
const char* s_dir_here      = "./";
const char* s_dot           = ".";
const char* s_symmetry      = "~";
const char* s_yes           = "yes";
const char* s_no            = "no";
const char* s_none          = "none";
const char* s_true          = ".true.";
const char* s_false         = ".false.";
const char* s_quote         = "\"";
const char* s_verbatim      = "\'";
const char* s_semicolon     = ";";
const char* s_colon         = ":";
const char* s_global        = "::";
const char* s_column        = ";";
const char* s_greater       = ">";
const char* s_less          = "<";
