#include "biddingga/initializers_2d.h"
#include <algorithm>
#include <vector>
#include "biddingga/helpers.h"

#include <eigen3/Eigen/Core>

using namespace Eigen;
using namespace numericaldists;
using namespace genericga;
using namespace genericga::binary;

namespace biddingga {

template <class Environment, class Phen>
std::vector<std::shared_ptr<multipop::SubGAAdapter<Environment, Phen>>>
MakeSubGAs(std::vector<Configuration2D> configs) {
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

SinglePopulationGA<ByteArrayGenotype, Grid> BinaryGA(
    int n_bits, const Configuration2D& config,
    std::function<Grid(const ByteArrayGenotype&)> phen_conv,
    std::function<std::vector<float>(const std::vector<Grid>&)> fit) {
  int n_bytes = (n_bits + CHAR_BIT - 1) / CHAR_BIT;
  auto genes = RandomGenes(config.n_strategies, n_bytes);
  Population<ByteArrayGenotype, Grid> init_pop(phen_conv, fit, genes);
  auto children_fact = std::make_unique<ChildrenFactory<ByteArrayGenotype>>(
      std::make_unique<SinglePointCrossover>(), std::make_unique<BitMutator>(2),
      std::make_unique<selector::TournamentMixed>(1));
  return SinglePopulationGA<ByteArrayGenotype, Grid>(
      std::move(init_pop), std::move(children_fact),
      std::make_unique<selector::ElitismDecorator>(
          std::make_unique<selector::TournamentMixed>(2.2), 2));
}

std::unique_ptr<CompositeGA<Grid>> MakeGAComposite(
    const Configuration2D& config) {
  assert(config.nx_segments % config.nx_composites == 0);
  assert(config.ny_segments % config.ny_composites == 0);
  int segs_per_xcomp = config.nx_segments / config.nx_composites;
  int n_floats_per_xcomp = segs_per_xcomp + 1;
  int segs_per_ycomp = config.ny_segments / config.ny_composites;
  int n_floats_per_ycomp = segs_per_ycomp + 1;
  int size = n_floats_per_xcomp * n_floats_per_ycomp;
  int n_bits = size * config.bit_precision;
  binary::FloatEncoding bid_enc(static_cast<float>(config.z_range.min),
                                static_cast<float>(config.z_range.max),
                                config.bit_precision);
  std::vector<binary::FloatEncoding> encodings(size, bid_enc);
  ArrayXd xvals = ArrayXd::LinSpaced(config.nx_segments + 1, config.x_range.min,
                                     config.x_range.max);
  ArrayXd yvals = ArrayXd::LinSpaced(config.ny_segments + 1, config.y_range.min,
                                     config.y_range.max);
  std::vector<std::shared_ptr<AbstractSinglePopulationGA<Grid>>> gas;
  for (int row = 0; row < config.ny_segments; row += segs_per_ycomp) {
    for (int col = 0; col < config.nx_segments; col += segs_per_xcomp) {
      ArrayXd sub_xvals = xvals.segment(col, n_floats_per_xcomp);
      ArrayXd sub_yvals = yvals.segment(row, n_floats_per_ycomp);
      std::function<Grid(const ByteArrayGenotype&)> conversion =
          GridSortDecorator{BinaryToGrid{sub_xvals, sub_yvals, encodings}};
      gas.push_back(
          std::make_shared<SinglePopulationGA<ByteArrayGenotype, Grid>>(
              BinaryGA(n_bits, config, conversion)));
    }
  }
  return std::make_unique<CompositeGA<Grid>>(
      gas, GridMeanJoiner{n_floats_per_xcomp, n_floats_per_ycomp,
                          config.nx_composites, config.ny_composites});
}

void SortByRow(ArrayXXd& arr) {
  for (int i = 0; i < arr.rows(); ++i) {
    ArrayXd row = arr.row(i);
    std::sort(row.data(), row.data() + row.size());
    arr.row(i) = row;
  }
}

Grid GridSortDecorator::operator()(const ByteArrayGenotype& gene) {
  Grid grid = phen_conv(gene);
  SortByRow(grid.zs);
  return grid;
}

Grid BinaryToGrid::operator()(const ByteArrayGenotype& gene) const {
  ArrayXf ys = gene.ToEigenFloatArray(nums_);
  ArrayXXf bids(1, ys.size());
  bids.row(0) = ys.segment(0, ys.size());
  bids.resize(valys_.size(), valxs_.size());
  return {valxs_, valys_, bids.cast<double>()};
}

PhenotypeStrategy2D GridJoinerSortDecorator::operator()(
    const std::vector<PhenotypeStrategy2D>& strats) {
  auto strat = joiner(strats);
  SortByRow(strat.phenotype.zs);
  return strat;
}

PhenotypeStrategy2D GridMeanJoiner::operator()(
    const std::vector<PhenotypeStrategy2D>& strats) {
  if (nx_chunks_ == 1 && ny_chunks_ == 1) {
    return strats[0];
  }
  assert(nx_chunks_ * ny_chunks_ == strats.size());
  PhenotypeStrategy2D joined{{JoinX(strats), JoinY(strats), JoinZ(strats)},
                             JoinFitnesses(strats)};
  return joined;
}

ArrayXd GridMeanJoiner::JoinX(const std::vector<PhenotypeStrategy2D>& strats) {
  int n_xs = nx_chunks_ * (xs_per_chunk_ - 1) + 1;
  ArrayXd joined_xs(n_xs);
  // Set endpoints
  joined_xs(0) = strats[0].phenotype.xs(0);
  joined_xs(n_xs - 1) = strats[nx_chunks_ - 1].phenotype.xs(xs_per_chunk_ - 1);
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

ArrayXd GridMeanJoiner::JoinY(const std::vector<PhenotypeStrategy2D>& strats) {
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

ArrayXXd GridMeanJoiner::JoinZ(const std::vector<PhenotypeStrategy2D>& strats) {
  int n_xs = nx_chunks_ * (xs_per_chunk_ - 1) + 1;
  int n_ys = ny_chunks_ * (ys_per_chunk_ - 1) + 1;
  ArrayXXd zs(n_ys, n_xs);

  // Set corners
  zs(0, 0) = strats[0].phenotype.zs(0, 0);
  zs(0, n_xs - 1) = strats[nx_chunks_ - 1].phenotype.zs(0, xs_per_chunk_ - 1);
  zs(n_ys - 1, 0) =
      strats[nx_chunks_ * (ny_chunks_ - 1)].phenotype.zs(ys_per_chunk_ - 1, 0);
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
          (strats[chunk_ind].phenotype.zs.col(0).segment(1, ys_per_chunk_ - 2) +
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
          (strats[chunk_ind].phenotype.zs.row(0).segment(1, xs_per_chunk_ - 2) +
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
      zs(j, i) =
          (strats[tr_ind].phenotype.zs(0, 0) +
           strats[tl_ind].phenotype.zs(0, xs_per_chunk_ - 1) +
           strats[br_ind].phenotype.zs(ys_per_chunk_ - 1, 0) +
           strats[bl_ind].phenotype.zs(ys_per_chunk_ - 1, xs_per_chunk_ - 1)) /
          4;
    }
  }

  return zs;
}

float GridMeanJoiner::JoinFitnesses(
    const std::vector<PhenotypeStrategy2D>& strats) {
  float fit = 0;
  for (const auto& strat : strats) {
    fit += strat.fitness;
  }
  return fit;
}
}  // namespace biddingga
