#ifndef GENERICGA_MULTIPOP_GA_H_
#define GENERICGA_MULTIPOP_GA_H_

#include <algorithm>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <utility>
#include <vector>
#include <iostream>

#include "genericga/multipop/abstract_sub_ga.h"

namespace genericga {
namespace multipop {

template <class Environment>
class GA {
 public:
  GA(std::vector<std::shared_ptr<AbstractSubGA<Environment>>> gas,
     Environment env, int batch_size = 1, int max_memory = 1)
      : gas_(std::move(gas)),
        gen_(std::random_device()()),
        batch_size_(batch_size),
        memory_(max_memory, env),
        max_memory_size_(max_memory),
        cur_memory_size_(0),
        mem_index_(0) {
    for (auto& ga : gas_) {
      ga->SubmitPlayStrat(memory_[mem_index_]);
    }
    IncrementMemory();
  }

  void RunRound(int n = 1) {
    for (int i = 0; i < n; ++i) {
      RunSingleRound();
    }
  }

 private:
  void IncrementMemory() {
    if (cur_memory_size_ < max_memory_size_) {
      ++cur_memory_size_;
    }
    mem_index_ = (mem_index_ + 1) % max_memory_size_;
  }
  void RunSingleRound() {
    std::vector<int> priorities(gas_.size());
    std::transform(gas_.begin(), gas_.end(), std::back_inserter(priorities),
                   [](const std::shared_ptr<AbstractSubGA<Environment>>& ga) {
                     return ga->GetPriority();
                   });
    std::set<int> unique_priorities(priorities.begin(), priorities.end());

    std::vector<Environment*> sel_envs;
    if (cur_memory_size_ < batch_size_) {
      for (int i = 0; i < cur_memory_size_; ++i) {
        sel_envs.push_back(&memory_[i]);
      }
    } else {
      std::vector<int> sel_inds(cur_memory_size_);
      std::iota(sel_inds.begin(), sel_inds.end(), 0);
      std::shuffle(sel_inds.begin(), sel_inds.end(), gen_);
      for (int i = 0; i < batch_size_; ++i) {
        sel_envs.push_back(&memory_[sel_inds[i]]);
      }
    }
    // All GAs of equal priority play strategies simulataneously.
    // Smaller priority moves first.  This can be useful when coordination is
    // required.
    for (auto p : unique_priorities) {
      // Update population.  (e.g. survey the environments and consider possible
      // moves)
      for (auto& ga : gas_) {
        if (ga->GetPriority() == p) {
          ga->RunRound(sel_envs);
        }
      }

      // Play strategy into each environment
      for (auto& ga : gas_) {
        if (ga->GetPriority() == p) {
          ga->SubmitPlayStrat(memory_[mem_index_]);
        }
      }
    }
    IncrementMemory();
  }

  std::mt19937 gen_;
  std::vector<std::shared_ptr<AbstractSubGA<Environment>>> gas_;
  std::vector<Environment> memory_;
  int batch_size_;
  int max_memory_size_;
  int cur_memory_size_ = 0;
  int mem_index_ = 0;
};

}  // namespace multipop
}  // namespace genericga

#endif  // GENERICGA_MULTIPOP_GA_H_
