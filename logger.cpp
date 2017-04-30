#ifndef logger_cpp
#define logger_cpp
#include <fstream>
#include <string>
#include <ctime>

using namespace std;

enum LogType { ERROR = -3, WARN = -2, OUTPUT = -1, INFO = 0, DEBUG = 1, TRACE = 2 };

class Logger {
public:
   LogType loglevel;
   ostream& of;
   bool echo;

  Logger() : of (cout) {
    //    of.open("log.txt");
    //    of = cout;
    loglevel = INFO;
    echo = false;
  }
  
   Logger( LogType inlevel, ostream& os, bool echo_in = false ): of( os ) {
    //of.open("log.txt");
    //    of = cout;
    loglevel = inlevel;
    echo = echo_in;
  }

  void set_level( LogType Lin ) {
    loglevel = Lin;
  }

  void operator()( LogType level, string msg ) {
    if (level <= loglevel) {
      switch (level) {
      case ERROR:
	of << time(NULL) << "\033[31m [ERROR] \033[0m" << msg << endl;
	break;
      case WARN:
	of << time(NULL) << "\033[33m [WARN] \033[0m" << msg << endl;
	break;
      case OUTPUT:
	of << time(NULL) << " [OUTPUT] " << msg << endl;
	break;
      case INFO:
	of << time(NULL) << "\033[32m [INFO] \033[0m" << msg << endl;
	break;
      case DEBUG:
	of << time(NULL) << "\033[31m [DEBUG] \033[0m" << msg << endl;
	break;
      case TRACE:
	of << time(NULL) << "\033[31m [TRACE] \033[0m" << msg << endl;
	break;
      }
      if (echo) {
	 if (level <= loglevel) {
	    switch (level) {
	    case ERROR:
	       cout << time(NULL) << "\033[31m [ERROR] \033[0m" << msg << endl;
	       break;
	    case WARN:
	       cout << time(NULL) << "\033[33m [WARN] \033[0m" << msg << endl;
	       break;
	    case OUTPUT:
	       cout << time(NULL) << " [OUTPUT] " << msg << endl;
	       break;
	    case INFO:
	       cout << time(NULL) << "\033[32m [INFO] \033[0m" << msg << endl;
	       break;
	    case DEBUG:
	       cout << time(NULL) << "\033[31m [DEBUG] \033[0m" << msg << endl;
	       break;
	    case TRACE:
	       cout << time(NULL) << "\033[31m [TRACE] \033[0m" << msg << endl;
	       break;
	    }
	 }
      }
    }	     
  }	     
  

  ~Logger() {
    //of.close();
  }
};

#endif
