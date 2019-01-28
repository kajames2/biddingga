#ifndef _NUMERICALDISTS_UNEVEN_BILERPER_
#define _NUMERICALDISTS_UNEVEN_BILERPER_

#include <functional>
#include <vector>

namespace numericaldists {

class UnevenBilerper {
 public:
  UnevenBilerper() {}
  UnevenBilerper(std::vector<float> xs, std::vector<float> ys,
                 std::vector<std::vector<float>> zs);
  UnevenBilerper(std::vector<float> ys, std::vector<std::function<float(float)>> slices);
  float operator()(float x, float y) const;

 private:
  int GetIndex(float y) const;
  float GetAlpha(int index, float y) const;
  std::vector<std::function<float(float)>> slices_;
  std::vector<float> ys_;
};

}  // namespace numericaldists

#endif  // _NUMERICALDISTS_UNEVEN_BILERPER_
