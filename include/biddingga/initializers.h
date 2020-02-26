#ifndef BIDDINGA_INITIALIZERS_H_
#define BIDDINGA_INITIALIZERS_H_

#include <algorithm>
#include <vector>
#include "configuration.h"

namespace biddingga {

std::unique_ptr<CompositeGA<numericaldists::Scatter>> MakeGAComposite(
    BidFunctionGAConfiguration config);

SinglePopulationGA<genericga::binary::ByteArrayGenotype,
                   numericaldists::Scatter>
BinaryGA(const Configuration1D& config,
         std::function<numericaldists::Scatter(
             const genericga::binary::ByteArrayGenotype&)>
             phen_conv,
         std::function<
             std::vector<float>(const std::vector<numericaldists::Scatter>&)>
             fit = [](const std::vector<numericaldists::Scatter>& phens) {
               return std::vector<float>(phens.size(), -1.0);
             });

struct BinaryToScatter {
  Scatter operator()(const genericga::binary::ByteArrayGenotype& gene) const;
  ArrayXd vals_;
  std::vector<binary::FloatEncoding> nums_;
};

struct ScatterSortDecorator {
  Scatter operator()(const genericga::binary::ByteArrayGenotype& gene) {
    Scatter scatter = phen_conv(gene);
    std::sort(scatter.ys.data(), scatter.ys.data() + scatter.ys.size());
    return scatter;
  }
  std::function<numericaldists::Scatter(
      const genericga::binary::ByteArrayGenotype& gene)>
      phen_conv;
};

struct ScatterReverseSortDecorator {
  Scatter operator()(const genericga::binary::ByteArrayGenotype& gene);
  std::function<numericaldists::Scatter(
      const genericga::binary::ByteArrayGenotype&)>
      phen_conv;
};

struct ScatterJoinerSortDecorator {
  PhenotypeStrategy<numericaldists::Scatter> operator()(
      const std::vector<PhenotypeStrategy<numericaldists::Scatter>>& strats);
  std::function<PhenotypeStrategy<numericaldists::Scatter>(
      const std::vector<PhenotypeStrategy<numericaldists::Scatter>>&)>
      joiner;
};

struct ScatterJoinerReverseSortDecorator {
  PhenotypeStrategy<numericaldists::Scatter> operator()(
      const std::vector<PhenotypeStrategy<numericaldists::Scatter>>& strats);
  std::function<PhenotypeStrategy<numericaldists::Scatter>(
      const std::vector<PhenotypeStrategy<numericaldists::Scatter>>&)>
      joiner;
};

struct ScatterSeparateJoiner {
  PhenotypeStrategy<numericaldists::Scatter> operator()(
      const std::vector<PhenotypeStrategy<numericaldists::Scatter>>& strats);
};

struct ScatterMedianJoiner {
  PhenotypeStrategy<numericaldists::Scatter> operator()(
      const std::vector<PhenotypeStrategy<numericaldists::Scatter>>& strats);
};
}  // namespace biddingga

#endif BIDDINGA_INITIALIZERS_H_
