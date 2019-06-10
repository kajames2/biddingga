#ifndef GENERICGA_PHENOTYPE_STRATEGY_H_
#define GENERICGA_PHENOTYPE_STRATEGY_H_

namespace genericga {

template <class Phen>
struct PhenotypeStrategy {
  Phen phenotype;
  float fitness = -1;
};

}  // namespace genericga
#endif  // GENERICGA_PHENOTYPE_STRATEGY_H_
