#ifndef NUMERICALDISTS_BOUNDS_H_
#define NUMERICALDISTS_BOUNDS_H_

#include <cstddef>
#include <ostream>

#include "numericaldists/interval.h"

namespace numericaldists {

struct Bounds {
  Bounds() {}
  Bounds(float a, float b) : a(a), b(b) {}
  explicit Bounds(Interval r) : a(r.min), b(r.max) {}
  float a;
  float b;
};

bool operator==(const Bounds& g1, const Bounds& g2);
std::ostream& operator<<(std::ostream& os, const Bounds& r);

}  // namespace numericaldists

#endif  // NUMERICALDISTS_BOUNDS_H_
