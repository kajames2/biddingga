#ifndef _GENERICGA_GENOTYPE_EVALUATOR_H_
#define _GENERICGA_GENOTYPE_EVALUATOR_H_

#include <algorithm>
#include <cassert>
#include <memory>
#include <vector>

#include "genericga/abstract_genotype_evaluator.h"
#include "genericga/phenotype_strategy.h"

namespace genericga {

template <class Gen, class Phen>
class GenotypeEvaluator : public AbstractGenotypeEvaluator<Gen> {
 public:
  GenotypeEvaluator(std::function<Phen(const Gen&)> phen_conv,
                    std::function<float(const Phen&)> fit_calc =
                        [](const Phen&) { return -1.0; });
  void SetFitnessCalculator(std::function<float(const Phen&)> fit_calc);
  float GetFitness(const Gen& gene) const override;
  std::vector<float> GetFitnesses(const std::vector<Gen>& genes) const override;
  Phen GetPhenotype(const Gen& gene) const;
  float GetPhenotypeFitness(const Phen& phen) const;
  PhenotypeStrategy<Phen> GetPhenotypeStrategy(const Gen& gene) const;

 private:
  std::function<Phen(const Gen&)> phen_conv_;
  std::function<float(const Phen&)> fit_calc_;
};

template <class Gen, class Phen>
GenotypeEvaluator<Gen, Phen>::GenotypeEvaluator(
    std::function<Phen(const Gen&)> phen_conv,
    std::function<float(const Phen&)> fit_calc)
    : phen_conv_(phen_conv), fit_calc_(fit_calc) {}

template <class Gen, class Phen>
void GenotypeEvaluator<Gen, Phen>::SetFitnessCalculator(
    std::function<float(const Phen&)> fit_calc) {
  fit_calc_ = fit_calc;
}

template <class Gen, class Phen>
Phen GenotypeEvaluator<Gen, Phen>::GetPhenotype(const Gen& gene) const {
  return phen_conv_(gene);
}

template <class Gen, class Phen>
float GenotypeEvaluator<Gen, Phen>::GetFitness(const Gen& gene) const {
  return GetPhenotypeFitness(GetPhenotype(gene));
}

template <class Gen, class Phen>
float GenotypeEvaluator<Gen, Phen>::GetPhenotypeFitness(
    const Phen& phen) const {
  return fit_calc_(phen);
}

template <class Gen, class Phen>
PhenotypeStrategy<Phen> GenotypeEvaluator<Gen, Phen>::GetPhenotypeStrategy(
    const Gen& gene) const {
  Phen phen = GetPhenotype(gene);
  return PhenotypeStrategy<Phen>{phen, GetPhenotypeFitness(phen)};
}

template <class Gen, class Phen>
std::vector<float> GenotypeEvaluator<Gen, Phen>::GetFitnesses(
    const std::vector<Gen>& genes) const {
  std::vector<float> fits;
  fits.reserve(genes.size());
  for (const Gen& gene : genes) {
    fits.emplace_back(GetFitness(gene));
  }
  return fits;
}

}  // namespace genericga
#endif  // _GENERICGA_GENOTYPE_EVALUATOR_H_
