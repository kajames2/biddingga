#include "numericaldists/bounds.h"

#include <cstddef>
#include <ostream>

#include <boost/container_hash/hash.hpp>

namespace numericaldists {

bool operator==(const Bounds& g1, const Bounds& g2) {
  return g1.a == g2.a && g1.b == g2.b;
}

std::size_t hash_value(const Bounds& s) {
  std::size_t seed = 0;
  boost::hash_combine(seed, s.a);
  boost::hash_combine(seed, s.b);

  return seed;
}

std::ostream& operator<<(std::ostream& os, const Bounds& r) {
  os << r.a << "->" << r.b;
  return os;
}

}  // namespace numericaldists
