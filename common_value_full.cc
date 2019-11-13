#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/uniform.hpp>

#include "auctions/common_value_signal_endpoints.h"

#include "genericga/binary/bit_mutator.h"
#include "genericga/binary/byte_array_genotype.h"
#include "genericga/binary/encoding.h"
#include "genericga/binary/single_point_crossover.h"
#include "genericga/composite_ga.h"
#include "genericga/multipop/ga.h"
#include "genericga/multipop/sub_ga_adapter.h"
#include "genericga/selector/elitism_decorator.h"
#include "genericga/selector/keep_best.h"
#include "genericga/selector/ranked_weighted.h"
#include "genericga/selector/tournament.h"
#include "genericga/selector/tournament_mixed.h"
#include "genericga/selector/tournament_poisson.h"
#include "genericga/single_population_ga.h"

#include "numericaldists/distribution.h"
#include "numericaldists/grid.h"

#include <algorithm>
#include <cstdlib>
#include <eigen3/Eigen/Core>
#include <iostream>
#include <vector>

using namespace genericga;
using namespace auctions;
using namespace numericaldists;
using namespace boost::math;
using namespace Eigen;

using Gen = binary::ByteArrayGenotype;
using Phen = Grid;

struct BidFunctionGAConfiguration {
  int id = -1;
  Interval valuex_range = {0, 1};
  Interval valuey_range = {0, 1};
  Interval bid_range = {0, 1};
  int nx_composites = 10;
  int ny_composites = 10;
  int n_strategies = 100;
  int n_children = 100;
  int nx_segments = 60;
  int ny_segments = 60;
  int bit_precision = 32;
};

template <class Phen>
SinglePopulationGA<binary::ByteArrayGenotype, Phen> BinaryGA(
    std::function<Phen(const binary::ByteArrayGenotype&)> phen_conv,
    int pop_size, int n_bits,
    std::function<std::vector<float>(const std::vector<Phen>&)> fit =
        [](const std::vector<Phen>& phens) {
          return std::vector<float>(phens.size(), -1.0);
        }) {
  using Gen = binary::ByteArrayGenotype;
  int n_bytes = (n_bits + CHAR_BIT - 1) / CHAR_BIT;
  std::vector<Gen> genes;
  auto generator = std::mt19937(std::random_device()());
  std::uniform_int_distribution<int> dist(0, UCHAR_MAX);
  for (int i = 0; i < pop_size; ++i) {
    std::vector<unsigned char> rand_gene(n_bytes);
    for (int j = 0; j < n_bytes; ++j) {
      rand_gene[j] = static_cast<unsigned char>(dist(generator));
    }
    genes.emplace_back(rand_gene);
  }
  Population<Gen, Phen> init_pop(phen_conv, fit, genes);

  auto children_fact = std::make_unique<ChildrenFactory<Gen>>(
      std::make_unique<binary::SinglePointCrossover>(),
      std::make_unique<binary::BitMutator>(2),
      std::make_unique<selector::TournamentMixed>(1));

  return SinglePopulationGA<binary::ByteArrayGenotype, Phen>(
      std::move(init_pop), std::move(children_fact),
      std::make_unique<selector::ElitismDecorator>(
          std::make_unique<selector::TournamentMixed>(2.2), 2));
}

void SortByRow(ArrayXXd& arr) {
  for (int i = 0; i < arr.rows(); ++i) {
    ArrayXd row = arr.row(i);
    std::sort(row.data(), row.data() + row.size());
    arr.row(i) = row;
  }
}

struct SortDecorator {
  Grid operator()(const Gen& gene) {
    Grid grid = phen_conv(gene);
    SortByRow(grid.zs);
    return grid;
  }
  std::function<Phen(const Gen& gene)> phen_conv;
};

struct BinaryToGrid {
  Grid operator()(const Gen& gene) const {
    ArrayXf ys = gene.ToEigenFloatArray(nums_);
    ArrayXXf bids(1, ys.size());
    bids.row(0) = ys.segment(0, ys.size());
    bids.resize(valys_.size(), valxs_.size());
    return {valxs_, valys_, bids.cast<double>()};
  }
  ArrayXd valxs_;
  ArrayXd valys_;
  std::vector<binary::FloatEncoding> nums_;
};

struct JoinerSortDecorator {
  PhenotypeStrategy<Phen> operator()(
      const std::vector<PhenotypeStrategy<Phen>>& strats) {
    auto strat = joiner(strats);
    SortByRow(strat.phenotype.zs);
    return strat;
  }
  std::function<PhenotypeStrategy<Phen>(
      const std::vector<PhenotypeStrategy<Phen>>&)>
      joiner;
};

struct MeanJoiner {
  PhenotypeStrategy<Grid> operator()(
      const std::vector<PhenotypeStrategy<Grid>>& strats) {
    if (nx_chunks_ == 1 && ny_chunks_ == 1) {
      return strats[0];
    }
    assert(nx_chunks_ * ny_chunks_ == strats.size());
    PhenotypeStrategy<Phen> joined{
        {JoinX(strats), JoinY(strats), JoinZ(strats)}, JoinFitnesses(strats)};
    return joined;
  }

  ArrayXd JoinX(const std::vector<PhenotypeStrategy<Grid>>& strats) {
    int n_xs = nx_chunks_ * (xs_per_chunk_ - 1) + 1;
    ArrayXd joined_xs(n_xs);
    // Set endpoints
    joined_xs(0) = strats[0].phenotype.xs(0);
    joined_xs(n_xs - 1) =
        strats[nx_chunks_ - 1].phenotype.xs(xs_per_chunk_ - 1);
    // Set interiors
    for (int i = 0; i < nx_chunks_; ++i) {
      joined_xs.segment(i * (xs_per_chunk_ - 1) + 1, xs_per_chunk_ - 2) =
          strats[i].phenotype.xs.segment(1, xs_per_chunk_ - 2);
    }
    // Set overlaps
    for (int i = 1; i < nx_chunks_; ++i) {
      joined_xs(i * (xs_per_chunk_ - 1)) =
          (strats[i - 1].phenotype.xs(xs_per_chunk_ - 1) +
           strats[i].phenotype.xs(0)) /
          2;
    }
    return joined_xs;
  }

  ArrayXd JoinY(const std::vector<PhenotypeStrategy<Grid>>& strats) {
    int n_ys = ny_chunks_ * (ys_per_chunk_ - 1) + 1;
    ArrayXd joined_ys(n_ys);
    // Set endpoints
    joined_ys(0) = strats[0].phenotype.ys(0);
    joined_ys(n_ys - 1) =
        strats[strats.size() - 1].phenotype.ys(ys_per_chunk_ - 1);
    // Set interiors
    for (int i = 0; i < ny_chunks_; ++i) {
      int chunk_ind = i * nx_chunks_;
      joined_ys.segment(i * (ys_per_chunk_ - 1) + 1, ys_per_chunk_ - 2) =
          strats[chunk_ind].phenotype.ys.segment(1, ys_per_chunk_ - 2);
    }
    // Set overlaps
    for (int i = 1; i < ny_chunks_; ++i) {
      int chunk_ind = i * nx_chunks_;
      joined_ys(i * (ys_per_chunk_ - 1)) =
          (strats[chunk_ind - 1].phenotype.ys(ys_per_chunk_ - 1) +
           strats[chunk_ind].phenotype.ys(0)) /
          2;
    }
    return joined_ys;
  }

  ArrayXXd JoinZ(const std::vector<PhenotypeStrategy<Grid>>& strats) {
    int n_xs = nx_chunks_ * (xs_per_chunk_ - 1) + 1;
    int n_ys = ny_chunks_ * (ys_per_chunk_ - 1) + 1;
    ArrayXXd zs(n_ys, n_xs);

    // Set corners
    zs(0, 0) = strats[0].phenotype.zs(0, 0);
    zs(0, n_xs - 1) = strats[nx_chunks_ - 1].phenotype.zs(0, xs_per_chunk_ - 1);
    zs(n_ys - 1, 0) = strats[nx_chunks_ * (ny_chunks_ - 1)].phenotype.zs(
        ys_per_chunk_ - 1, 0);
    zs(n_ys - 1, n_xs - 1) = strats[strats.size() - 1].phenotype.zs(
        ys_per_chunk_ - 1, xs_per_chunk_ - 1);

    // Set top/bottom interiors
    for (int i = 0; i < nx_chunks_; ++i) {
      int chunk_ind = i;
      zs.row(0).segment(i * (xs_per_chunk_ - 1) + 1, xs_per_chunk_ - 2) =
          strats[chunk_ind].phenotype.zs.row(0).segment(1, xs_per_chunk_ - 2);
      chunk_ind += nx_chunks_ * (ny_chunks_ - 1);
      zs.row(n_ys - 1).segment(i * (xs_per_chunk_ - 1) + 1, xs_per_chunk_ - 2) =
          strats[chunk_ind]
              .phenotype.zs.row(ys_per_chunk_ - 1)
              .segment(1, xs_per_chunk_ - 2);
    }

    // Set right/left interiors
    for (int i = 0; i < ny_chunks_; ++i) {
      int chunk_ind = i * nx_chunks_;
      zs.col(0).segment(i * (ys_per_chunk_ - 1) + 1, ys_per_chunk_ - 2) =
          strats[chunk_ind].phenotype.zs.col(0).segment(1, ys_per_chunk_ - 2);
      chunk_ind += nx_chunks_ - 1;
      zs.col(n_xs - 1).segment(i * (ys_per_chunk_ - 1) + 1, ys_per_chunk_ - 2) =
          strats[chunk_ind]
              .phenotype.zs.col(xs_per_chunk_ - 1)
              .segment(1, ys_per_chunk_ - 2);
    }

    // Set top/bottom midpoints
    for (int i = 1; i < nx_chunks_; ++i) {
      int chunk_ind = i;
      int overlap_ind = chunk_ind - 1;
      zs(0, i * (xs_per_chunk_ - 1)) =
          (strats[chunk_ind].phenotype.zs(0, 0) +
           strats[overlap_ind].phenotype.zs(0, xs_per_chunk_ - 1)) /
          2;
      chunk_ind += nx_chunks_ * (ny_chunks_ - 1);
      overlap_ind = chunk_ind - 1;
      zs(n_ys - 1, i * (xs_per_chunk_ - 1)) =
          (strats[chunk_ind].phenotype.zs(ys_per_chunk_ - 1, 0) +
           strats[overlap_ind].phenotype.zs(ys_per_chunk_ - 1,
                                            xs_per_chunk_ - 1)) /
          2;
    }

    // Set left/right boundary midpoints
    for (int i = 1; i < ny_chunks_; ++i) {
      int chunk_ind = i * nx_chunks_;
      int overlap_ind = chunk_ind - nx_chunks_;
      zs(i * (ys_per_chunk_ - 1), 0) =
          (strats[chunk_ind].phenotype.zs(0, 0) +
           strats[overlap_ind].phenotype.zs(ys_per_chunk_ - 1, 0)) /
          2;
      chunk_ind += nx_chunks_ - 1;
      overlap_ind = chunk_ind - nx_chunks_;
      zs(i * (ys_per_chunk_ - 1), n_xs - 1) =
          (strats[chunk_ind].phenotype.zs(0, xs_per_chunk_ - 1) +
           strats[overlap_ind].phenotype.zs(ys_per_chunk_ - 1,
                                            xs_per_chunk_ - 1)) /
          2;
    }

    // Set true interiors
    for (int row = 0; row < ny_chunks_; ++row) {
      for (int col = 0; col < nx_chunks_; ++col) {
        int chunk_ind = row * nx_chunks_ + col;
        int i = col * (xs_per_chunk_ - 1) + 1;
        int j = row * (ys_per_chunk_ - 1) + 1;
        zs.block(j, i, ys_per_chunk_ - 2, xs_per_chunk_ - 2) =
            strats[chunk_ind].phenotype.zs.block(1, 1, ys_per_chunk_ - 2,
                                                 xs_per_chunk_ - 2);
      }
    }

    // Set interior x-overlaps
    for (int row = 0; row < ny_chunks_; ++row) {
      for (int col = 1; col < nx_chunks_; ++col) {
        int chunk_ind = row * nx_chunks_ + col;
        int overlap_ind = chunk_ind - 1;
        int i = col * (xs_per_chunk_ - 1);
        int j = row * (ys_per_chunk_ - 1) + 1;
        zs.col(i).segment(j, ys_per_chunk_ - 2) =
            (strats[chunk_ind].phenotype.zs.col(0).segment(1,
                                                           ys_per_chunk_ - 2) +
             strats[overlap_ind]
                 .phenotype.zs.col(xs_per_chunk_ - 1)
                 .segment(1, ys_per_chunk_ - 2)) /
            2;
      }
    }

    // Set interior y-overlaps
    for (int row = 1; row < ny_chunks_; ++row) {
      for (int col = 0; col < nx_chunks_; ++col) {
        int chunk_ind = row * nx_chunks_ + col;
        int overlap_ind = chunk_ind - nx_chunks_;
        int i = col * (xs_per_chunk_ - 1) + 1;
        int j = row * (ys_per_chunk_ - 1);
        zs.row(j).segment(i, xs_per_chunk_ - 2) =
            (strats[chunk_ind].phenotype.zs.row(0).segment(1,
                                                           xs_per_chunk_ - 2) +
             strats[overlap_ind]
                 .phenotype.zs.row(ys_per_chunk_ - 1)
                 .segment(1, xs_per_chunk_ - 2)) /
            2;
      }
    }

    // Set interior x&y/overlaps
    for (int row = 1; row < ny_chunks_; ++row) {
      for (int col = 1; col < nx_chunks_; ++col) {
        int tr_ind = row * nx_chunks_ + col;
        int tl_ind = tr_ind - 1;
        int br_ind = tr_ind - nx_chunks_;
        int bl_ind = tr_ind - nx_chunks_ - 1;
        int i = col * (xs_per_chunk_ - 1);
        int j = row * (ys_per_chunk_ - 1);
        zs(j, i) = (strats[tr_ind].phenotype.zs(0, 0) +
                    strats[tl_ind].phenotype.zs(0, xs_per_chunk_ - 1) +
                    strats[br_ind].phenotype.zs(ys_per_chunk_ - 1, 0) +
                    strats[bl_ind].phenotype.zs(ys_per_chunk_ - 1,
                                                xs_per_chunk_ - 1)) /
                   4;
      }
    }

    return zs;
  }

  float JoinFitnesses(const std::vector<PhenotypeStrategy<Grid>>& strats) {
    float fit = 0;
    for (const auto& strat : strats) {
      fit += strat.fitness;
    }
    return fit;
  }
  int xs_per_chunk_;
  int ys_per_chunk_;
  int nx_chunks_;
  int ny_chunks_;
};

template <class Phen>
std::unique_ptr<CompositeGA<Phen>> MakeGAComposite(
    BidFunctionGAConfiguration config) {
  assert(config.nx_segments % config.nx_composites == 0);
  assert(config.ny_segments % config.ny_composites == 0);
  int segs_per_xcomp = config.nx_segments / config.nx_composites;
  int n_floats_per_xcomp = segs_per_xcomp + 1;
  int segs_per_ycomp = config.ny_segments / config.ny_composites;
  int n_floats_per_ycomp = segs_per_ycomp + 1;
  int size = 2 * (n_floats_per_xcomp * n_floats_per_ycomp);
  int n_bits = size * config.bit_precision;
  binary::FloatEncoding bid_enc(static_cast<float>(config.bid_range.min),
                                static_cast<float>(config.bid_range.max),
                                config.bit_precision);
  std::vector<binary::FloatEncoding> encodings(size, bid_enc);
  ArrayXd xvals = ArrayXd::LinSpaced(
      config.nx_segments + 1, config.valuex_range.min, config.valuex_range.max);
  ArrayXd yvals = ArrayXd::LinSpaced(
      config.ny_segments + 1, config.valuey_range.min, config.valuey_range.max);
  std::vector<std::shared_ptr<AbstractSinglePopulationGA<Phen>>> gas;
  for (int row = 0; row < config.ny_segments; row += segs_per_ycomp) {
    for (int col = 0; col < config.nx_segments; col += segs_per_xcomp) {
      ArrayXd sub_xvals = xvals.segment(col, n_floats_per_xcomp);
      ArrayXd sub_yvals = yvals.segment(row, n_floats_per_ycomp);
      std::function<Phen(const Gen&)> conversion =
          SortDecorator{BinaryToGrid{sub_xvals, sub_yvals, encodings}};
      gas.push_back(std::make_shared<SinglePopulationGA<Gen, Phen>>(
          BinaryGA<Phen>(conversion, config.n_strategies, n_bits)));
    }
  }
  return std::make_unique<CompositeGA<Phen>>(
      gas, MeanJoiner{n_floats_per_xcomp, n_floats_per_ycomp,
                      config.nx_composites, config.ny_composites});
}

template <class Environment, class Phen>
std::vector<std::shared_ptr<multipop::SubGAAdapter<Environment, Phen>>>
MakeSubGAs(std::vector<BidFunctionGAConfiguration> configs) {
  std::vector<std::shared_ptr<multipop::SubGAAdapter<Environment, Phen>>>
      ind_gas;
  for (int i = 0; i < configs.size(); ++i) {
    auto ga = MakeGAComposite<Phen>(std::move(configs[i]));
    auto subga = std::make_shared<multipop::SubGAAdapter<Environment, Phen>>(
        std::move(ga), i, 0, std::make_unique<selector::KeepBest>());
    ind_gas.push_back(subga);
  }
  return ind_gas;
}

template <class Environment, class Phen>
multipop::GA<Environment> MakeMultipopDriver(
    const std::vector<
        std::shared_ptr<multipop::SubGAAdapter<Environment, Phen>>>& gas,
    Environment env) {
  std::vector<std::shared_ptr<multipop::AbstractSubGA<Environment>>> sub_gas(
      gas.begin(), gas.end());
  return multipop::GA<Environment>(sub_gas, env);
}

template <class Environment, class Phen>
void RunAndOutput(
    multipop::GA<Environment>& driver,
    std::vector<std::shared_ptr<multipop::SubGAAdapter<Environment, Phen>>>&
        gas,
    Environment& env, int n_rounds, int output_frequency) {
  for (int n = 0; n < n_rounds / output_frequency; ++n) {
    driver.RunRound(output_frequency);
    for (int j = 0; j < gas.size(); ++j) {
      std::cout << gas[j]->GetBestStrategy().phenotype << std::endl;
      env.AcceptStrategy(gas[j]->GetBestStrategy().phenotype, j);
    }
    for (int j = 0; j < gas.size(); ++j) {
      std::cout << env.GetFitness(gas[j]->GetBestStrategy().phenotype, j);
      if (j != gas.size() - 1) {
        std::cout << ",";
      }
    }
    std::cout << std::endl;
  }
}

int main(int argc, char** argv) {
  std::vector<int> n_draws;
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " Draws1, Draws2 ..." << std::endl;
    return 1;
  } else {
    for (int i = 1; i < argc; ++i) {
      n_draws.push_back(std::atoi(argv[i]));
    }
  }

  Distribution value_dist = uniform_distribution<>(500, 9500);
  Distribution error_dist = uniform_distribution<>(-500, 500);
  std::vector<BidFunctionGAConfiguration> configs;
  for (int i = 0; i < n_draws.size(); ++i) {
    BidFunctionGAConfiguration config;
    config.id = i;
    config.valuex_range = {lower(value_dist) + lower(error_dist),
                           upper(value_dist) + upper(error_dist)};
    config.valuey_range = {0, 2 * (upper(error_dist) - lower(error_dist))};
    config.bid_range = {lower(value_dist) + lower(error_dist),
                        upper(value_dist) + upper(error_dist)};
    configs.push_back(std::move(config));
  }

  CommonValueSignalEndpoints auction(value_dist, error_dist, n_draws);
  auto gas = MakeSubGAs<CommonValueSignalEndpoints, Grid>(configs);
  auto driver =
      MakeMultipopDriver<CommonValueSignalEndpoints, Grid>(gas, auction);

  int n_rounds = 3000;
  int output_frequency = 30;
  RunAndOutput(driver, gas, auction, n_rounds, output_frequency);

  return 0;
}
