#include "numericaldists/line_segment.h"

#include <cstddef>
#include <ostream>

#include <boost/container_hash/hash.hpp>

#include "numericaldists/bounds.h"
#include "numericaldists/interval.h"

namespace numericaldists {

LineSegment::LineSegment(Interval xint, Bounds ybounds)
    : xint_(xint),
      ybounds_(ybounds),
      inv_interval_(1 / (xint.max - xint.min)),
      slope_((ybounds.a - ybounds.b) / (xint.max - xint.min)) {}

LineSegment::LineSegment(Bounds xbounds, Bounds ybounds) {
  if (xbounds.a <= xbounds.b) {
    xint_ = Interval{xbounds.a, xbounds.b};
    ybounds_ = ybounds;
  } else {
    xint_ = Interval(xbounds.b, xbounds.a);
    ybounds_ = Bounds(ybounds.b, ybounds.a);
  }
  inv_interval_ = 1 / (xint_.max - xint_.min);
  slope_ = (ybounds_.a - ybounds_.b) / (xint_.max - xint_.min);
}

float LineSegment::GetY(float x) const {
  float alpha = (x - xint_.min) * inv_interval_;
  return alpha * ybounds_.b + (1 - alpha) * ybounds_.a;
}

bool operator==(const LineSegment& g1, const LineSegment& g2) {
  return g1.xint_ == g2.xint_ && g1.ybounds_ == g2.ybounds_;
}

std::size_t hash_value(const LineSegment& s) {
  std::size_t seed = 0;
  boost::hash_combine(seed, s.xint_);
  boost::hash_combine(seed, s.ybounds_);

  return seed;
}

std::ostream& operator<<(std::ostream& os, const LineSegment& ls) {
  os << ls.GetXInterval() << "\t:\t" << ls.GetYBounds();
  return os;
}

}  // namespace numericaldists
