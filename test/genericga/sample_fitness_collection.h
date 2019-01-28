#ifndef _SAMPLE_FITNESS_COLLECTION_H_
#define _SAMPLE_FITNESS_COLLECTION_H_

#include <vector>

#include "genericga/fitness_collection.h"

namespace gatests {

class SampleFitnessCollection : public genericga::FitnessCollection {
 public:
  SampleFitnessCollection(std::vector<float> vec = std::vector<float>{1, 3, 6,
                                                                      -2})
      : fits(vec) {}
  std::vector<float> GetFitnesses() const override { return fits; }

 protected:
  float GetFitness(int i) const override { return fits[i]; }

 private:
  std::vector<float> fits;
};

}  // namespace gatests

#endif  // _SAMPLE_FITNESS_COLLECTION_H_
