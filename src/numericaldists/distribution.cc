#include "numericaldists/distribution.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace numericaldists {

float lower(const Distribution& dist) {
  float q = quantile(dist, 0);
  return std::isinf(q) ? quantile(dist, 0.000001) : q;
}
float upper(const Distribution& dist) {
  float q = quantile(dist, 1);
  return std::isinf(q) ? quantile(dist, 0.999999) : q;
}

float lower(const std::vector<Distribution>& dists) {
  return lower(*std::min_element(
      dists.begin(), dists.end(),
      [](const Distribution& dist1, const Distribution& dist2) {
        return lower(dist1) < lower(dist2);
      }));
}

float upper(const std::vector<Distribution>& dists) {
  return upper(*std::max_element(
      dists.begin(), dists.end(),
      [](const Distribution& dist1, const Distribution& dist2) {
        return upper(dist1) < upper(dist2);
      }));
}

}  // namespace numericaldists
