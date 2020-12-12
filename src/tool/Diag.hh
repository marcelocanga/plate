#ifndef __DIAG_MEC_HH__
#define __DIAG_MEC_HH__
//#pragma interface
#include <cstdio>
#include <cstdarg>
#include <iostream>
#include <fstream>
#include "Essential.hh"
//
//************************************************      external
//
class Diag;
class Ostring;
class Analyzer;
class String;

extern Diag diag;
//
//************************************************      
//
struct DiagPrecision{
  int p;
};

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//
//               -----  class Diag:public ostream  -----
//
//					    Copyright
//                                          Marcelo Canga
//					    Year 2000-2010
//                                          All Rights Reserved
//
//
// C: none: no messages
//    help: ask for help on commands 
//    fatal: error with program termination. May or may not be a programming error. Usually is.
//    program: logic programming error. should be fatal but we let go. requires developer assistance
//    error: error in computations or choices, invalidates run
//    warning: not an error but a warning
//    caution: a warning that cannot be cleared, that is not registered like one
//    info: important information messages,
//    echo: print information during debugging, should not be permanent
//    debug: information for debugging purposes, debug2,debug3 same but more details
//    main: important modules,
//    inout: important  info inside subroutines,
//    detail: details inside main subroutines,
//    function: details inside auxiliary functions
//    factory:  error at the factory level
//    variable: parser variable def, i.e., name,word,integer,etc,
//    type: messages inside string, array, etc.
//
//    For error, exception, fatal error and warnings use
//    diag_mesg(diag_error,<<vars)  (not constants)
//    I check whether the number of error exceded the max
// 
//    For debug use debug. This is the latest, and meant to be transformed or deleted
//    For &command -help use help. This should always be emmited
//
//    To use the log facitility for things that obvious need printing
//    use help
//
//    For debug,output,etc,use
//    if(diag.level(diag.inout))
//        diag<<var1<<etc.
//    It is type safe, it is fast, it is cool but I cannot control
//    what to print after the command is executed.
//    Even better, the macro,
//    message(Diag::warning,a<<b<<c<<d<<etc); will in addition add
//    the top and bottom messages for error, warnings
//    
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

class Diag{

public:

  enum  DiagLevel  {security,none,help,fatal,program,error,fail,warning,caution,info,echo,debug,debug1,debug2,debug3,main,inout,detail,
                    function,factory,variable,type};
  enum DiagLevel   diag_error_compatible;

  enum MosesCompatible { moses_diaglevel_o = 01, moses_inertia_o = 02, moses_emit_o = 04, moses_language_o = 010 };  // set with &debug

  enum ExitCode  {normal_e,warning_e,fail_e,error_e,fatal_e,program_e,security_e,layout_e,restarted_e};

  enum ThrowCode {fail_t,security_t,fatal_t,program_t,max_error_t,
		  layout_t,                                         // layout execution mode;
		  none_t} throw_code;
  
  bool  x_echo, r_echo, m_echo;
  bool  x_echo_encr, r_echo_encr, m_echo_encr;

private:
  enum DiagLevel loglevel,termlevel, diaglevel, severity, form_level;
  bool is_report_status;            // no, means do not report errors, useful for tests that have error, yet finish normal
  bool is_encripted;
  bool is_restarted;

  int   nerror,nwarn,nfail,ncaution, max_error;


public:

  friend class Command;

  std::ofstream fos;
  std::ostream  &stdos;
  std::ostream  &stderror;
  std::stringstream stdss;

  int moses_compatible_flag;
//
//************************************************      creation,del
//
  Diag();
  Diag(char* name);
  void nmax_error(int nmax);
  void nmax_warn(int nmax);
  int get_nerror(){ return nerror; }
  int get_nwarn() { return nwarn; }
  int get_ncaution(){ return ncaution; }
  int get_nfail() { return nfail; }
  void  initialize();
  ~Diag();
  std::filebuf* rdbuf();
//
//*********************************************      security
//
  void set_encripted(bool val);
  bool get_encripted(){return is_encripted;}
//
//************************************************      device operations
//
  void open(const char* nm);
  void close();
//
//************************************************      output operations
//
  void mesg (enum DiagLevel, const char*, ...);
  void message (enum DiagLevel, const char* ,va_list&);
  void send_message (enum DiagLevel, const char* ,va_list&);
  void send_message (enum DiagLevel, const char* , ...);
  void top_message (enum DiagLevel);
  void bot_message (enum DiagLevel);
  void status();
  int  exit_code();
  void set_code(enum ThrowCode);

  void boot_interface();
  void diag_interface(Analyzer& an);
  void common_block_size_interface();
  bool moses_flag(enum MosesCompatible moses_flag){ return moses_flag & moses_compatible_flag; }
  bool spl_flag(int);
  void send_diag(std::stringstream&);
//
//************************************************      device state
//
  void setloglevel(enum DiagLevel dl);
  void settermlevel(enum DiagLevel dl);
  void setlevel(enum DiagLevel dl);
  
  void set_mecho(bool pd);
  void set_recho(bool pd);
  void set_xecho(bool pd);
  void set_restarted(bool val);

  bool test_level(enum DiagLevel sseverity);
  bool term_level(enum DiagLevel sseverity);
  
//
//************************************************      form
//
  void form( std::string fmt, ...);
  void oform( std::string fmt, ...);
  void vform(std::ostream& os, std::string fmt, va_list& ap);
//
//*********************************************      broadcast
//
  void broadcast(String& msg);
  void broadcast(double time,double dtime, double stop);
//
//************************************************      operators
//
  template< class T>

  Diag& operator<<(T& p){

    stdss.str(std::string());
    
    if(severity <= termlevel){
      if(severity > Diag::info){
	stderror << p;
      }
      else{
	stdss    << p;
	send_diag(stdss);
      }
    }
    
    if(severity <= loglevel )
        fos<<p;

    return *this; 
  }
//
//*********************************************      precision setter
//
  DiagPrecision& precision(int prec)
  {
    static DiagPrecision os;
    os.p = prec;
    return os;
  } 
//
//*********************************************      
//
  Diag& operator<<(Ostring& p);              // use for formated string
  Diag& operator<<(const char* p);
  Diag& operator<<(const char p);
  Diag& operator<<(DiagPrecision);
  Diag& operator<<(std::ostream& omanip(std::ostream &));
  Diag& operator<<(String& );
//
//************************************************      done
//
};

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//
//               -----  inline void Diag::mesg  -----
//
//					    Copyright
//                                          Marcelo Canga
//					    Year 2000-2004
//                                          All Rights Reserved
//
//
// C: inline message call
//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

inline void Diag::mesg(enum DiagLevel severity, const char *fmt, ...)
{
//
//************************************************      Quick check
//
  if(!test_level(severity))
      return;
//
//************************************************      parse message
//
  va_list ap;
  va_start(ap,fmt);
  message(severity,fmt,ap);
  va_end(ap);
//
//************************************************      done
//
}


#endif
