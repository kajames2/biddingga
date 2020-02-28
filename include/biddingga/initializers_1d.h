#ifndef BIDDINGA_INITIALIZERS_1D_H_
#define BIDDINGA_INITIALIZERS_1D_H_

#include "genericga/binary/byte_array_genotype.h"
#include "genericga/composite_ga.h"
#include "genericga/multipop/ga.h"
#include "genericga/multipop/sub_ga_adapter.h"
#include "genericga/phenotype_strategy.h"
#include "genericga/single_population_ga.h"
#include "numericaldists/scatter.h"

#include <eigen3/Eigen/Core>

#include <vector>
#include "configuration.h"

namespace biddingga {

using PhenotypeStrategy1D =
    genericga::PhenotypeStrategy<numericaldists::Scatter>;

std::unique_ptr<genericga::CompositeGA<numericaldists::Scatter>>
MakeGAComposite(const Configuration1D& config);

genericga::SinglePopulationGA<genericga::binary::ByteArrayGenotype,
                              numericaldists::Scatter>
BinaryGA(int n_bits, const Configuration1D& config,
         std::function<numericaldists::Scatter(
             const genericga::binary::ByteArrayGenotype&)>
             phen_conv,
         std::function<
             std::vector<float>(const std::vector<numericaldists::Scatter>&)>
             fit = [](const std::vector<numericaldists::Scatter>& phens) {
               return std::vector<float>(phens.size(), -1.0);
             });

struct BinaryToScatter {
  numericaldists::Scatter operator()(
      const genericga::binary::ByteArrayGenotype& gene) const;
  Eigen::ArrayXd vals_;
  std::vector<genericga::binary::FloatEncoding> nums_;
};

struct ScatterSortDecorator {
  numericaldists::Scatter operator()(
      const genericga::binary::ByteArrayGenotype& gene);
  std::function<numericaldists::Scatter(
      const genericga::binary::ByteArrayGenotype& gene)>
      phen_conv;
};

struct ScatterReverseSortDecorator {
  numericaldists::Scatter operator()(
      const genericga::binary::ByteArrayGenotype& gene);
  std::function<numericaldists::Scatter(
      const genericga::binary::ByteArrayGenotype&)>
      phen_conv;
};

struct ScatterJoinerSortDecorator {
  PhenotypeStrategy1D operator()(
      const std::vector<PhenotypeStrategy1D>& strats);
  std::function<PhenotypeStrategy1D(const std::vector<PhenotypeStrategy1D>&)>
      joiner;
};

struct ScatterJoinerReverseSortDecorator {
  PhenotypeStrategy1D operator()(
      const std::vector<PhenotypeStrategy1D>& strats);
  std::function<PhenotypeStrategy1D(const std::vector<PhenotypeStrategy1D>&)>
      joiner;
};

struct ScatterSeparateJoiner {
  PhenotypeStrategy1D operator()(
      const std::vector<PhenotypeStrategy1D>& strats);
};

struct ScatterMedianJoiner {
  PhenotypeStrategy1D operator()(
      const std::vector<PhenotypeStrategy1D>& strats);
};
}  // namespace biddingga

#endif  // BIDDINGA_INITIALIZERS_1D_H_
