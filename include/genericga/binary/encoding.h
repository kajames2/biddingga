#ifndef _GENERICGA_BINARY_ENCODING_H_
#define _GENERICGA_BINARY_ENCODING_H_

namespace genericga {
namespace binary {

struct Encoding {
  int bit_precision;
  float min;
  float max;
  bool is_gray_coded = true;
};

}  // namespace binary
}  // namespace genericga

#endif  // _GENERICGA_BINARY_ENCODING_H_
