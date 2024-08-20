#include <ginac/ginac.h>
#include <stdexcept>

// Note: we use just simple power series implementations of Li functions.
// These could be replaced by simple Fortran code instead of linking to GiNaC.
// However, using GiNaC has the advantage to be able to switch to high
// precision evaluations if needed.

using std::runtime_error;
using GiNaC::ex;
using GiNaC::numeric;
using GiNaC::sqrt;
using GiNaC::Li;
using GiNaC::log;
using GiNaC::is_a;
using GiNaC::ex_to;
using GiNaC::zeta;
using GiNaC::lst;
using GiNaC::Pi;

#ifdef DIGITS

typedef ex numtype;
namespace {
static inline ex nn(int a) { return numeric(a); }
static inline ex nn(int a, int b) { return numeric(a,b); }
static inline ex tonn(const ex& x) { return x.evalf(); }
static inline double todouble(const ex& x) {
  ex n = x.evalf();
  if (!is_a<numeric>(n))  throw runtime_error("can't evaluate polylog");
  return ex_to<numeric>(n).to_double();
}
}

#else

typedef double numtype;
namespace {
static inline double nn(long long a) { return static_cast<double>(a); }
static inline double nn(long long a, long long b) { return a/static_cast<double>(b); }
static inline double tonn(const ex& x) {
  ex n = x.evalf();
  if (!is_a<numeric>(n))  throw runtime_error("can't evaluate polylog");
  return ex_to<numeric>(n).to_double();
}
static inline double todouble(const ex& x) { return tonn(x); }
}

#endif

