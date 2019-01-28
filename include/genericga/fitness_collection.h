#ifndef _GENERICGA_FITNESS_COLLECTION_H_
#define _GENERICGA_FITNESS_COLLECTION_H_

#include "genericga/selector.h"

#include <vector>
#include <cassert>

// A collection of fitness values with tools for ranking while preserving the
// underlying collection order.
namespace genericga {
class FitnessCollection {
 public:
  virtual std::vector<float> GetFitnesses() const = 0;
  // Ranking fitnesses.  Ties get the average of the ranks at that value.
  // ex. fitnesses of 5,6,5,1,5 have ranks: 2,4,2,0,2 (2 = mean(1,2,3))
  std::vector<float> GetFitnessRankings() const;

  // First element is index of smallest element, next is index of 2nd smallest
  // element, etc. ex. fitnesses of 5,8,1 have orderings: 2,0,1
  std::vector<int> GetFitnessOrderings() const;
  virtual int Size() const;

 protected:
  virtual float GetFitness(int i) const = 0;
  int SelectIndex(Selector& selector) const {
    return SelectIndices(selector, 1)[0];
  }
  
  std::vector<int> SelectIndices(Selector& selector, int n) const {
    assert(Size() > 0 && "Cannot select from empty collection.");
    return selector.SelectIndices(*this, n);
  }
};

}  // namespace genericga

#endif  // _GENERICGA_FITNESS_COLLECTION_H_
