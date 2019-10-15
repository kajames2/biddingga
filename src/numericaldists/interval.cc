#include "numericaldists/interval.h"

#include <cstddef>
#include <ostream>

namespace numericaldists {

bool InInterval(Interval interval, double val) {
  return val >= interval.min && val <= interval.max;
}

double GetSpan(Interval interval) { return interval.max - interval.min; }

bool operator==(const Interval& g1, const Interval& g2) {
  return g1.min == g2.min && g1.max == g2.max;
}

std::ostream& operator<<(std::ostream& os, const Interval& r) {
  os << '[' << r.min << ',' << r.max << ')';
  return os;
}

}  // namespace numericaldists
