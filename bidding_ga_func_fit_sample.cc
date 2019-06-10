#include "genericga/binary/bit_mutator.h"
#include "genericga/binary/byte_array_genotype.h"
#include "genericga/binary/encoding.h"
#include "genericga/binary/single_point_crossover.h"
#include "genericga/selector/elitism_decorator.h"
#include "genericga/selector/keep_best.h"
#include "genericga/selector/ranked_weighted.h"
#include "genericga/selector/tournament.h"
#include "genericga/selector/tournament_mixed.h"
#include "genericga/selector/tournament_poisson.h"
#include "genericga/single_population_ga.h"

#include <omp.h>
#include <cmath>

using namespace genericga;

using Gen = binary::ByteArrayGenotype;
using Phen = std::vector<float>;

struct BinaryGAConfiguration {
  std::unique_ptr<Crossover<binary::ByteArrayGenotype>> crossover =
      std::make_unique<binary::SinglePointCrossover>();
  std::unique_ptr<Mutator<binary::ByteArrayGenotype>> mutation =
      std::make_unique<binary::BitMutator>(2);
  std::unique_ptr<Selector> parent_selector =
      std::make_unique<selector::TournamentMixed>(1);
  std::unique_ptr<Selector> survival_selector =
      std::make_unique<selector::ElitismDecorator>(
          std::make_unique<selector::TournamentMixed>(2.2), 2);
  int n_strategies = 500;
  int n_children = 500;
};

std::vector<binary::ByteArrayGenotype> RandomBinaryStrategies(int pop_size,
                                                              int n_bits) {
  int n_bytes = (n_bits + CHAR_BIT - 1) / CHAR_BIT;
  std::vector<binary::ByteArrayGenotype> genes;
  auto generator = std::mt19937(std::random_device()());
  std::uniform_int_distribution<int> dist(0, UCHAR_MAX);
  for (int i = 0; i < pop_size; ++i) {
    std::vector<unsigned char> rand_gene(n_bytes);
    for (int j = 0; j < n_bytes; ++j) {
      rand_gene[j] = static_cast<unsigned char>(dist(generator));
    }
    genes.emplace_back(rand_gene);
  }
  return genes;
}

SinglePopulationGA<binary::ByteArrayGenotype, std::vector<float>> BinaryGA(
    std::vector<binary::FloatEncoding> encodings,
    std::function<float(const std::vector<float>&)> fit,
    BinaryGAConfiguration config = BinaryGAConfiguration()) {
  int n_bits = 0;
  for (const auto& e : encodings) {
    n_bits += e.bit_precision;
  }
  auto init_genes = RandomBinaryStrategies(config.n_strategies, n_bits);
  auto phen_conv =
      [encodings](const binary::ByteArrayGenotype& gene) -> std::vector<float> {
    return gene.ToFloatArray(encodings);
  };

  auto par_fit = [fit](const std::vector<std::vector<float>>& strats) {
    std::vector<float> fits(strats.size());
#pragma omp parallel for
    for (int i = 0; i < strats.size(); ++i) {
      fits[i] = fit(strats[i]);
    }
    return fits;
  };

  Population<binary::ByteArrayGenotype, std::vector<float>> init_pop(
      phen_conv, par_fit, init_genes);

  auto children_fact =
      std::make_unique<ChildrenFactory<binary::ByteArrayGenotype>>(
          std::move(config.crossover), std::move(config.mutation),
          std::move(config.parent_selector));

  return SinglePopulationGA<binary::ByteArrayGenotype, std::vector<float>>(
      std::move(init_pop), std::move(children_fact),
      std::move(config.survival_selector));
}

float SineAccuracy(std::vector<float> coeffs) {
  float sqr_err = 0;
  for (float x = -3; x < 3; x += 0.01) {
    float est = 0;
    for (int i = 0; i < coeffs.size(); ++i) {
      est += coeffs[i] * std::pow(x, i);
    }
    sqr_err += std::pow(est - std::sin(x), 2);
  }
  return -sqr_err;
}

int main() {
  auto ga = BinaryGA({4, binary::FloatEncoding{-20, 20}}, SineAccuracy);
  ga.RunRound(1000);
  auto coeffs = ga.GetBestStrategy().phenotype;
  for (auto c : coeffs) {
    std::cout << c << ",";
  }
  std::cout << std::endl;
  return 0;
}
