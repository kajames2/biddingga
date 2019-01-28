#ifndef _NUMERICALDISTS_LINE_SEGMENT_H_
#define _NUMERICALDISTS_LINE_SEGMENT_H_

#include <cstddef>
#include <iostream>

#include "numericaldists/bounds.h"
#include "numericaldists/interval.h"

namespace numericaldists {

class LineSegment {
 public:
  LineSegment(Interval xint, Bounds yint);
  LineSegment(Bounds xbounds, Bounds ybounds);
  float GetY(float x) const;
  Interval GetXInterval() const { return xint_; }
  Bounds GetYBounds() const { return ybounds_; }
  LineSegment GetInverseLineSegment() const {
    return LineSegment(ybounds_, Bounds(xint_));
  }
  float GetSlope() const { return slope_; }
  bool IsInInterval(float x) const { return InInterval(xint_, x); }
  friend std::size_t hash_value(const LineSegment& s);
  friend bool operator==(const LineSegment& g1, const LineSegment& g2);

 private:
  Interval xint_;
  Bounds ybounds_;
  float inv_interval_;
  float slope_;
};

bool operator==(const LineSegment& g1, const LineSegment& g2);
std::size_t hash_value(const LineSegment& s);
std::ostream& operator<<(std::ostream& os, const LineSegment& ls);

}  // namespace numericaldists

#endif  // _NUMERICALDISTS_LINE_SEGMENT_H_
