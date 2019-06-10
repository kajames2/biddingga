#ifndef GENERICGA_MUTATOR_H_
#define GENERICGA_MUTATOR_H_

namespace genericga {

template <class Gen>
class Mutator {
 public:
  virtual void operator()(Gen& genotype) = 0;
  virtual ~Mutator() {}
};
}  // namespace genericga

#endif  // GENERICGA_MUTATOR_H_
