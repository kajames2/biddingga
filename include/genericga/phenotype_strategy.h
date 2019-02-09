#ifndef _GENERICGA_PHENOTYPE_STRATEGY_H_
#define _GENERICGA_PHENOTYPE_STRATEGY_H_

namespace genericga {

template <class Phen>
struct PhenotypeStrategy {
  Phen phenotype;
  float fitness = -1;
};

}  // namespace genericga
#endif  // _GENERICGA_PHENOTYPE_STRATEGY_H_
