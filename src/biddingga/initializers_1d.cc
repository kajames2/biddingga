#include "biddingga/initializers_1d.h"
#include <algorithm>
#include <vector>
#include "biddingga/helpers.h"

#include <eigen3/Eigen/Core>

using namespace genericga::binary;
using namespace numericaldists;
using namespace Eigen;
using namespace genericga;

namespace biddingga {

template <class Environment, class Phen>
std::vector<std::shared_ptr<multipop::SubGAAdapter<Environment, Phen>>>
MakeSubGAs(std::vector<Configuration1D> configs) {
  std::vector<std::shared_ptr<multipop::SubGAAdapter<Environment, Phen>>>
      ind_gas;
  for (int i = 0; i < configs.size(); ++i) {
    auto ga = MakeGAComposite(std::move(configs[i]));
    auto subga = std::make_shared<multipop::SubGAAdapter<Environment, Phen>>(
        std::move(ga), i, 0, std::make_unique<selector::KeepBest>());
    ind_gas.push_back(subga);
  }
  return ind_gas;
}

SinglePopulationGA<ByteArrayGenotype, Scatter> BinaryGA(
    int n_bits, const Configuration1D& config,
    std::function<Scatter(const ByteArrayGenotype&)> phen_conv,
    std::function<std::vector<float>(const std::vector<Scatter>&)> fit) {
  int n_bytes = (n_bits + CHAR_BIT - 1) / CHAR_BIT;
  auto genes = RandomGenes(config.n_strategies, n_bytes);
  Population<ByteArrayGenotype, Scatter> init_pop(phen_conv, fit, genes);
  auto children_fact = std::make_unique<ChildrenFactory<ByteArrayGenotype>>(
      std::make_unique<SinglePointCrossover>(), std::make_unique<BitMutator>(2),
      std::make_unique<selector::TournamentMixed>(1));
  return SinglePopulationGA<ByteArrayGenotype, Scatter>(
      std::move(init_pop), std::move(children_fact),
      std::make_unique<selector::ElitismDecorator>(
          std::make_unique<selector::TournamentMixed>(2.2), 2));
}

std::unique_ptr<CompositeGA<Scatter>> MakeGAComposite(
    const Configuration1D& config) {
  assert(config.nx_segments % config.nx_composites == 0);
  int segs_per_comp = config.nx_segments / config.nx_composites;
  int n_floats_per_comp = segs_per_comp + 1;
  int n_bits = n_floats_per_comp * config.bit_precision;
  FloatEncoding bid_enc(static_cast<float>(config.y_range.min),
                        static_cast<float>(config.y_range.max),
                        config.bit_precision);
  std::vector<FloatEncoding> encodings(n_floats_per_comp, bid_enc);
  ArrayXd vals = ArrayXd::LinSpaced(config.nx_segments + 1, config.x_range.min,
                                    config.x_range.max);
  std::vector<std::shared_ptr<AbstractSinglePopulationGA<Scatter>>> gas;
  for (int i = 0; i < config.nx_segments; i += segs_per_comp) {
    ArrayXd sub_vals = vals.segment(i, n_floats_per_comp);
    std::function<Scatter(const ByteArrayGenotype&)> conversion =
        BinaryToScatter{sub_vals, encodings};
    gas.push_back(
        std::make_shared<SinglePopulationGA<ByteArrayGenotype, Scatter>>(
            BinaryGA(n_bits, config, conversion)));
  }
  return std::make_unique<CompositeGA<Scatter>>(gas, ScatterMedianJoiner());
}

Scatter ScatterSortDecorator::operator()(const ByteArrayGenotype& gene) {
  Scatter scatter = phen_conv(gene);
  std::sort(scatter.ys.data(), scatter.ys.data() + scatter.ys.size());
  return scatter;
}

Scatter ScatterReverseSortDecorator::operator()(const ByteArrayGenotype& gene) {
  Scatter scatter = phen_conv(gene);
  std::sort(scatter.ys.data(), scatter.ys.data() + scatter.ys.size());
  std::reverse(scatter.ys.data(), scatter.ys.data() + scatter.ys.size());
  return scatter;
}

Scatter BinaryToScatter::operator()(const ByteArrayGenotype& gene) const {
  ArrayXf ys = gene.ToEigenFloatArray(nums_);
  return {vals_, ys.cast<double>()};
}

PhenotypeStrategy1D ScatterJoinerSortDecorator::operator()(
    const std::vector<PhenotypeStrategy1D>& strats) {
  auto strat = joiner(strats);
  std::sort(strat.phenotype.ys.data(),
            strat.phenotype.ys.data() + strat.phenotype.ys.size());
  return strat;
}

PhenotypeStrategy1D ScatterJoinerReverseSortDecorator::operator()(
    const std::vector<PhenotypeStrategy1D>& strats) {
  auto strat = joiner(strats);
  std::sort(strat.phenotype.ys.data(),
            strat.phenotype.ys.data() + strat.phenotype.ys.size());
  std::reverse(strat.phenotype.ys.data(),
               strat.phenotype.ys.data() + strat.phenotype.ys.size());
  return strat;
}

PhenotypeStrategy1D ScatterSeparateJoiner::operator()(
    const std::vector<PhenotypeStrategy1D>& strats) {
  int n_strats = strats.size();
  int strat_size = strats[0].phenotype.xs.size();
  ArrayXd joined_xs(strat_size * n_strats);
  ArrayXd joined_ys(strat_size * n_strats);
  float fit = 0;
  for (int i = 0; i < n_strats; ++i) {
    fit += strats[i].fitness;
    joined_xs.segment(strat_size * i, strat_size) << strats[i].phenotype.xs;
    joined_ys.segment(strat_size * i, strat_size) << strats[i].phenotype.ys;
  }
  return PhenotypeStrategy1D{{joined_xs, joined_ys}, fit};
}

PhenotypeStrategy1D ScatterMedianJoiner::operator()(
    const std::vector<PhenotypeStrategy1D>& strats) {
  int n_strats = strats.size();
  int strat_size = strats[0].phenotype.xs.size();
  ArrayXd joined_xs(n_strats * (strat_size - 1) + 1);
  ArrayXd joined_ys(n_strats * (strat_size - 1) + 1);
  float fit = 0;
  joined_xs(0) = strats[0].phenotype.xs(0);
  joined_ys(0) = strats[0].phenotype.ys(0);
  for (int i = 0; i < n_strats - 1; ++i) {
    fit += strats[i].fitness;
    joined_xs.segment(i * (strat_size - 1) + 1, strat_size - 1)
        << strats[i].phenotype.xs.segment(1, strat_size - 1);
    joined_ys.segment(i * (strat_size - 1) + 1, strat_size - 2)
        << strats[i].phenotype.ys.segment(1, strat_size - 2);
    joined_ys((i + 1) * (strat_size - 1)) =
        (strats[i].phenotype.ys(strat_size - 1) +
         strats[i + 1].phenotype.ys(0)) /
        2;
  }
  fit += strats[n_strats - 1].fitness;
  joined_xs.segment((n_strats - 1) * (strat_size - 1) + 1, strat_size - 1)
      << strats[strats.size() - 1].phenotype.xs.segment(1, strat_size - 1);
  joined_ys.segment((n_strats - 1) * (strat_size - 1) + 1, strat_size - 1)
      << strats[strats.size() - 1].phenotype.ys.segment(1, strat_size - 1);
  return PhenotypeStrategy1D{{joined_xs, joined_ys}, fit};
}

}  // namespace biddingga
