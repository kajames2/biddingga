#ifndef GENERICGA_REAL_SINGLE_POINT_CROSSOVER_H_
#define GENERICGA_REAL_SINGLE_POINT_CROSSOVER_H_

#include <climits>
#include <iostream>
#include <random>
#include <vector>

#include "genericga/crossover.h"

namespace genericga {
namespace real {

class SinglePointCrossover : public Crossover<std::vector<double>> {
 public:
  SinglePointCrossover() : gen_(std::random_device()()), dist(0, 0) {}

  void operator()(std::vector<double>& gene1,
                  std::vector<double>& gene2) override {
    if (gene1.size() < gene_size_ || gene2.size() < gene_size_) {
      gene_size_ = std::min(gene1.size(), gene2.size());
      dist = std::uniform_int_distribution<>(0, gene_size_ - 1);
    }

    int ind = dist(gen_);
    for (int i = ind; i < gene_size_; ++i) {
      std::swap(gene1[i], gene2[i]);
    }
  }

 private:
  std::mt19937 gen_;
  int gene_size_ = 0;
  std::uniform_int_distribution<> dist;
};

}  // namespace real
}  // namespace genericga

#endif  // GENERICGA_REAL_SINGLE_POINT_CROSSOVER_H_
