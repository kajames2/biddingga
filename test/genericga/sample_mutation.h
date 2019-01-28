#ifndef _SAMPLE_MUTATION_H_
#define _SAMPLE_MUTATION_H_

#include "genericga/mutator.h"

namespace gatests {

class SampleMutation
    : public genericga::Mutator<int> {
public:
  void operator()(int &gen) override {
    gen = gen * gen - 1;
  }
};

} // namespace gatests

#endif // _SAMPLE_MUTATION_H_


