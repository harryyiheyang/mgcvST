#include <Rcpp.h>
#include <algorithm>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef PSAPI_VERSION
#define PSAPI_VERSION 2
#endif
#include <windows.h>
#include <psapi.h>
#endif

// Read counters only; never allocate memory to test availability.
// [[Rcpp::export]]
Rcpp::NumericVector mgcvst_memory_status_cpp() {
  double available = NA_REAL, resident = NA_REAL;
#ifdef _WIN32
  MEMORYSTATUSEX status;
  status.dwLength = sizeof(status);
  if (GlobalMemoryStatusEx(&status)) {
    available = static_cast<double>(status.ullAvailPhys);
  }
  PROCESS_MEMORY_COUNTERS counters;
  if (GetProcessMemoryInfo(GetCurrentProcess(), &counters, sizeof(counters))) {
    resident = static_cast<double>(counters.WorkingSetSize);
  }
#endif
  return Rcpp::NumericVector::create(
    Rcpp::Named("available") = available, Rcpp::Named("resident") = resident);
}
