#ifndef _NUMERICALDISTS_UNEVEN_PIECEWISE_LINEAR_
#define _NUMERICALDISTS_UNEVEN_PIECEWISE_LINEAR_

#include <ostream>
#include <vector>

#include "numericaldists/interval.h"
#include "numericaldists/line_segment.h"

namespace numericaldists {

class UnevenPiecewiseLinear {
 public:
  UnevenPiecewiseLinear() {}
  virtual ~UnevenPiecewiseLinear() {}
  UnevenPiecewiseLinear(const std::vector<float>& xs, const std::vector<float>& ys);
  virtual LineSegment GetLine(float x) const;
  float GetBid(float x) const;
  std::vector<float> GetBids(std::vector<float> xs) const;
  friend std::size_t hash_value(const UnevenPiecewiseLinear& s);
  friend bool operator==(const UnevenPiecewiseLinear& g1,
                         const UnevenPiecewiseLinear& g2);
  friend std::ostream& operator<<(std::ostream& os,
                                  const UnevenPiecewiseLinear& pl);
  UnevenPiecewiseLinear(std::vector<LineSegment> segments);
  float operator()(float x) const { return GetBid(x); }
  std::vector<float> operator()(std::vector<float> xs) const { return GetBids(xs); }
  std::vector<LineSegment> GetSegments() const { return segments_; }
 protected:
  std::vector<LineSegment> segments_;
  Interval interval_;
};
bool operator==(const UnevenPiecewiseLinear& g1,
                const UnevenPiecewiseLinear& g2);
std::size_t hash_value(const UnevenPiecewiseLinear& s);

}  // namespace numericaldists

#endif  // _NUMERICALDISTS_UNEVEN_PIECEWISE_LINEAR_
