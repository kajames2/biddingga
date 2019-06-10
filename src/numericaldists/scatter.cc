#include "numericaldists/scatter.h"

namespace numericaldists {

std::ostream& operator<<(std::ostream& os, const Scatter& points) {
  for (int i = 0; i < points.xs.size() - 1; ++i) {
    os << points.xs(i) << ",";
  }
  os << points.xs(points.xs.size() -1) << "\n";
  for (int i = 0; i < points.ys.size() - 1; ++i) {
    os << points.ys(i) << ",";
  }
  os << points.ys(points.ys.size() -1);
  return os;
}

}  // namespace numericaldists
