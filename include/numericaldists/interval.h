#ifndef _NUMERICALDISTS_INTERVAL_H_
#define _NUMERICALDISTS_INTERVAL_H_

#include <cassert>
#include <cstddef>
#include <ostream>

namespace numericaldists {

struct Interval {
  Interval() {}
  Interval(float min, float max) : min(min), max(max) { assert(min <= max); }
  float min;
  float max;
};

bool InInterval(Interval interval, float val);
float GetSpan(Interval interval);
bool operator==(const Interval& g1, const Interval& g2);
std::size_t hash_value(const Interval& s);
std::ostream& operator<<(std::ostream& os, const Interval& r);

}  // namespace numericaldists

#endif  // _NUMERICALDISTS_INTERVAL_H_
