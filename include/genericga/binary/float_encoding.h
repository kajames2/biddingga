#ifndef GENERICGA_BINARY_FLOAT_ENCODING_H_
#define GENERICGA_BINARY_FLOAT_ENCODING_H_

namespace genericga {
namespace binary {

struct FloatEncoding {
  FloatEncoding(float min, float max, int bit_precision = 32,
                bool gray_coded = true)
      : min(min),
        max(max),
        bit_precision(bit_precision),
        is_gray_coded(gray_coded) {}
  float min;
  float max;
  int bit_precision = 32;
  bool is_gray_coded = true;
};

}  // namespace binary
}  // namespace genericga

#endif  // GENERICGA_BINARY_FLOAT_ENCODING_H_
