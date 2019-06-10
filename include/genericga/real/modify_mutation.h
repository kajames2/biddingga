#ifndef GENERICGA_REAL_MODIFY_MUTATION_H_
#define GENERICGA_REAL_MODIFY_MUTATION_H_

#include <random>
#include <vector>

#include "genericga/mutator.h"

namespace genericga {
namespace real {

class ModifyMutation : public Mutator<std::vector<double>> {
 public:
  ModifyMutation(double prob_rate, double std_dev, double min, double max)
      : gen_(std::random_device()()),
        muts_dist_(prob_rate),
        ind_dist_(0, 0),
        val_dist_(0, std_dev),
        min_(min),
        max_(max),
        gene_size_(0) {}

  void operator()(std::vector<double>& gene) override {
    if (gene.size() != gene_size_) {
      gene_size_ = gene.size();
      ind_dist_ = std::uniform_int_distribution<>(0, gene_size_ - 1);
    }

    int n_muts = muts_dist_(gen_);
    for (int i = 0; i < n_muts; ++i) {
      int ind = ind_dist_(gen_);
      gene[ind] = Mutate(gene[ind]);
    }
  }

 private:
  double Mutate(double val) {
    double new_val = val + val_dist_(gen_);
    if (new_val < min_ || new_val > max_) {
      return Mutate(val);
    }
    return new_val;
  }
  std::mt19937 gen_;
  std::poisson_distribution<> muts_dist_;
  std::uniform_int_distribution<> ind_dist_;
  std::normal_distribution<> val_dist_;
  int gene_size_;
  double min_;
  double max_;
};
}  // namespace real
}  // namespace genericga

#endif  // GENERICGA_REAL_MODIFY_MUTATION_H_
