#ifndef _GENERICGA_PHENOTYPE_STRATEGY_H_
#define _GENERICGA_PHENOTYPE_STRATEGY_H_

#include <boost/functional/hash.hpp>

namespace genericga {

template <class Phen>
struct PhenotypeStrategy {
  Phen phenotype;
  float fitness = -1;
};

template <class Phen>
bool operator<(const PhenotypeStrategy<Phen>& s1,
               const PhenotypeStrategy<Phen>& s2) {
  return s1.phenotype < s2.phenotype;
}

}  // namespace genericga
#endif  // _GENERICGA_PHENOTYPE_STRATEGY_H_
