#ifndef BIDDINGA_INITIALIZERS_2D_H_
#define BIDDINGA_INITIALIZERS_2D_H_

#include "biddingga/configuration.h"
#include "genericga/binary/byte_array_genotype.h"
#include "genericga/composite_ga.h"
#include "genericga/multipop/ga.h"
#include "genericga/multipop/sub_ga_adapter.h"
#include "genericga/phenotype_strategy.h"
#include "genericga/single_population_ga.h"
#include "numericaldists/grid.h"
#include "numericaldists/scatter.h"

#include <eigen3/Eigen/Core>

#include <vector>
#include "configuration.h"

namespace biddingga {

using PhenotypeStrategy2D = genericga::PhenotypeStrategy<numericaldists::Grid>;

std::unique_ptr<genericga::CompositeGA<numericaldists::Grid>> MakeGAComposite(
    const Configuration2D& config);

genericga::SinglePopulationGA<genericga::binary::ByteArrayGenotype,
                              numericaldists::Grid>
BinaryGA(
    int n_bits, const Configuration2D& config,
    std::function<
        numericaldists::Grid(const genericga::binary::ByteArrayGenotype&)>
        phen_conv,
    std::function<std::vector<float>(const std::vector<numericaldists::Grid>&)>
        fit = [](const std::vector<numericaldists::Grid>& phens) {
          return std::vector<float>(phens.size(), -1.0);
        });

struct GridSortDecorator {
  numericaldists::Grid operator()(
      const genericga::binary::ByteArrayGenotype& gene);
  std::function<numericaldists::Grid(
      const genericga::binary::ByteArrayGenotype& gene)>
      phen_conv;
};

struct BinaryToGrid {
  numericaldists::Grid operator()(
      const genericga::binary::ByteArrayGenotype& gene) const;
  Eigen::ArrayXd valxs_;
  Eigen::ArrayXd valys_;
  std::vector<genericga::binary::FloatEncoding> nums_;
};

struct GridJoinerSortDecorator {
  PhenotypeStrategy2D operator()(
      const std::vector<PhenotypeStrategy2D>& strats);
  std::function<PhenotypeStrategy2D(const std::vector<PhenotypeStrategy2D>&)>
      joiner;
};

struct GridMeanJoiner {
  PhenotypeStrategy2D operator()(
      const std::vector<PhenotypeStrategy2D>& strats);
  Eigen::ArrayXd JoinX(const std::vector<PhenotypeStrategy2D>& strats);
  Eigen::ArrayXd JoinY(const std::vector<PhenotypeStrategy2D>& strats);
  Eigen::ArrayXXd JoinZ(const std::vector<PhenotypeStrategy2D>& strats);
  float JoinFitnesses(const std::vector<PhenotypeStrategy2D>& strats);
  int xs_per_chunk_;
  int ys_per_chunk_;
  int nx_chunks_;
  int ny_chunks_;
};

}  // namespace biddingga

#endif  // BIDDINGA_INITIALIZERS_2D_H_
