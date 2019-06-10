#ifndef GENERICGA_CROSSOVER_H_
#define GENERICGA_CROSSOVER_H_

namespace genericga {

//  Mixes two genotypes together.
template <class Gen>
class Crossover {
 public:
  virtual void operator()(Gen& genotype1, Gen& genotype2) = 0;
  virtual ~Crossover() {}
};
}  // namespace genericga

#endif  // GENERICGA_CROSSOVER_H_
