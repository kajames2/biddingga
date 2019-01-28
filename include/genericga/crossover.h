#ifndef _GENERICGA_CROSSOVER_H_
#define _GENERICGA_CROSSOVER_H_

namespace genericga {

//  Mixes two genotypes together.
template <class Gen>
class Crossover {
 public:
  virtual void operator()(Gen& genotype1, Gen& genotype2);
  virtual ~Crossover() {}
};
}  // namespace genericga

#endif  // _GENERICGA_CROSSOVER_H_
