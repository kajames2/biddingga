#ifndef _GENERICGA_MUTATOR_H_
#define _GENERICGA_MUTATOR_H_

namespace genericga {

template <class Gen>
class Mutator {
 public:
  virtual void operator()(Gen &genotype);
  virtual ~Mutator() {}
};
}  // namespace genericga

#endif  // _GENERICGA_MUTATOR_H_
