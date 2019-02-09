#ifndef _SAMPLE_FITNESS_CALCULATOR_H_
#define _SAMPLE_FITNESS_CALCULATOR_H_

#include <vector>

namespace gatests {

class SampleFitnessCalculator {
 public:
  explicit SampleFitnessCalculator(int n) : n_(n) {}
  std::vector<float> operator()(const std::vector<int>& phens_in) const {
    std::vector<float> phens(phens_in.begin(), phens_in.end());
    for (float& phen : phens) {
      phen += n_;
    }
    return phens;
  }

 private:
  int n_;
};

}  // namespace gatests

#endif  // _SAMPLE_FITNESS_CALCULATOR_H_
