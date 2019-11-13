#include "numericaldists/grid.h"

#include <ostream>

namespace numericaldists {

std::ostream& operator<<(std::ostream& os, const Grid& grid) {
  for (int i = 0; i < grid.xs.size() - 1; ++i) {
    os << grid.xs(i) << ",";
  }
  os << grid.xs(grid.xs.size() - 1) << "\n";
  for (int i = 0; i < grid.ys.size() - 1; ++i) {
    os << grid.ys(i) << ",";
  }
  os << grid.ys(grid.ys.size() - 1) << "\n";

  for (int row = 0; row < grid.ys.size(); ++row) {
    for (int col = 0; col < grid.xs.size(); ++col) {
      os << grid.zs(row, col);
      if (row != grid.ys.size() - 1 || col != grid.xs.size() - 1) {
        os << ",";
      }
    }
  }
  return os;
}

}  // namespace numericaldists
