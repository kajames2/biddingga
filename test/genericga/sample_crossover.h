#ifndef _SAMPLE_CROSSOVER_H_
#define _SAMPLE_CROSSOVER_H_

#include "genericga/crossover.h"

namespace gatests {

class SampleCrossover : public genericga::Crossover<int> {
 public:
  void operator()(int &gen1, int &gen2) override {
    int avg = (gen1 + gen2) / 2;
    gen1 = avg;
    gen2 = avg;
  }
};

}  // namespace gatests

#endif  // _SAMPLE_CROSSOVER_H_
