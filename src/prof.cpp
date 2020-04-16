#include "prof.h"
#ifdef DO_PROF
#include <gperftools/profiler.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <sstream>
#include <string>
#include <Rcpp.h>
#endif

#ifdef DO_PROF
std::atomic<bool> profiler::running_profiler(false);
#endif

profiler::profiler(const std::string &name)
{
#ifdef DO_PROF
  if(running_profiler)
    Rcpp::stop("Already running profiler...");
  running_profiler = true;

  std::stringstream ss;
  using namespace std::chrono;
  auto now = system_clock::now();
  auto ms = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;
  auto timer = system_clock::to_time_t(now);
  std::tm bt = *std::localtime(&timer);

  ss << "profile-" << name << std::put_time(&bt, "-%d-%m-%Y-%H-%M-%S-")
     << std::setfill('0') << std::setw(3) << ms.count()
     << ".log";
  Rcpp::Rcout << "Saving profile output to '" << ss.str() << "'" << std::endl;
  const std::string s = ss.str();
  ProfilerStart(s.c_str());
#endif
}

#ifdef DO_PROF
profiler::~profiler(){
  ProfilerStop();
  running_profiler = false;
}
#else
profiler::~profiler() = default;
#endif
