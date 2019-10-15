#ifndef NUMERICALDISTS_INTERVAL_H_
#define NUMERICALDISTS_INTERVAL_H_

#include <cassert>
#include <cstddef>
#include <ostream>

namespace numericaldists {

struct Interval {
  Interval() {}
  Interval(double min, double max) : min(min), max(max) { assert(min <= max); }
  double min;
  double max;
};

bool InInterval(Interval interval, double val);
double GetSpan(Interval interval);
bool operator==(const Interval& g1, const Interval& g2);
std::ostream& operator<<(std::ostream& os, const Interval& r);

}  // namespace numericaldists

#endif  // NUMERICALDISTS_INTERVAL_H_
