#include "numericaldists/interval.h"

#include <cstddef>
#include <ostream>

#include <boost/container_hash/hash.hpp>

namespace numericaldists {

bool InInterval(Interval interval, float val) {
  return val >= interval.min && val < interval.max;
}

float GetSpan(Interval interval) { return interval.max - interval.min; }

bool operator==(const Interval& g1, const Interval& g2) {
  return g1.min == g2.min && g1.max == g2.max;
}

std::size_t hash_value(const Interval& s) {
  std::size_t seed = 0;
  boost::hash_combine(seed, s.min);
  boost::hash_combine(seed, s.max);

  return seed;
}

std::ostream& operator<<(std::ostream& os, const Interval& r) {
  os << '[' << r.min << ',' << r.max << ')';
  return os;
}

}  // namespace numericaldists
