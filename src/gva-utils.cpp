#include "gva-utils.h"

namespace GaussHermite {
namespace GVA {

double mlogit_integral
  (double const mu, double const sigma, HermiteData<double> const &hd){
  return mlogit_integral_atomic<double>::comp(mu, M_SQRT2 * sigma, hd);
}

} // namespace GaussHermite
} // namespace GVA
