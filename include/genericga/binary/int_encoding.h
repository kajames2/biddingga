#ifndef GENERICGA_BINARY_INT_ENCODING_H_
#define GENERICGA_BINARY_INT_ENCODING_H_

#include <cmath>

namespace genericga {
namespace binary {

struct IntEncoding {
  IntEncoding(int min, int max, bool gray_coded = true)
      : min(min),
        max(max),
        is_gray_coded(gray_coded),
        n_bits(static_cast<int>(std::log2(max - min + 1) + 1)) {}
  int min;
  int max;
  int n_bits;
  bool is_gray_coded = true;
};

}  // namespace binary
}  // namespace genericga

#endif  // GENERICGA_BINARY_INT_ENCODING_H_
