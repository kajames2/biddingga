#ifndef _SAMPLE_FITNESS_CALCULATOR_H_
#define _SAMPLE_FITNESS_CALCULATOR_H_

namespace gatests {

class SampleFitnessCalculator {
public:
  explicit SampleFitnessCalculator(int n) : n_(n) {}
  float operator()(const int &phenotype) const {
    return phenotype + n_;
  }

private:
  int n_;
};

} // namespace gatests

#endif // _SAMPLE_FITNESS_CALCULATOR_H_
