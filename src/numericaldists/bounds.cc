#include "numericaldists/bounds.h"

#include <cstddef>
#include <ostream>

namespace numericaldists {

bool operator==(const Bounds& g1, const Bounds& g2) {
  return g1.a == g2.a && g1.b == g2.b;
}

std::ostream& operator<<(std::ostream& os, const Bounds& r) {
  os << r.a << "->" << r.b;
  return os;
}

}  // namespace numericaldists
