#ifndef _GENERICGA_MULTIPOP_GA_H_
#define _GENERICGA_MULTIPOP_GA_H_

#include <memory>
#include <set>
#include <vector>
#include "genericga/multipop/abstract_sub_ga.h"

namespace genericga {
namespace multipop {

template <class Environment>
class GA {
 public:
  GA(std::vector<std::shared_ptr<AbstractSubGA<Environment>>> gas,
     Environment env, int n_selections = 5)
      : gas_(std::move(gas)),
        envs_(n_selections, env),
        n_selections_(n_selections) {
    for (auto& ga : gas_) {
      for (auto& env : envs_) {
        ga->SubmitPlayStrat(env);
      }
    }
  }

  void RunRound(int n = 1) {
    for (int i = 0; i < n; ++i) {
      RunSingleRound();
    }
  }

 private:
  void RunSingleRound() {
    std::vector<int> priorities(gas_.size());
    std::transform(gas_.begin(), gas_.end(), std::back_inserter(priorities),
                   [](const std::shared_ptr<AbstractSubGA<Environment>>& ga) {
                     return ga->GetPriority();
                   });
    std::set<int> unique_priorities(priorities.begin(), priorities.end());

    // All GAs of equal priority play strategies simulataneously.
    // Smaller priority moves first.  This can be useful when coordination is
    // required.
    for (auto p : unique_priorities) {
      // Update population.  (e.g. survey the environments and consider possible
      // moves)
      for (auto& ga : gas_) {
        if (ga->GetPriority() == p) {
          ga->RunRound(envs_);
        }
      }

      // Play strategy into the environments
      for (auto& ga : gas_) {
        if (ga->GetPriority() == p) {
          for (auto& env : envs_) {
            ga->SubmitPlayStrat(env);
          }
        }
      }
    }
  }

  std::vector<std::shared_ptr<AbstractSubGA<Environment>>> gas_;
  std::vector<Environment> envs_;
  int n_selections_;
};

}  // namespace multipop
}  // namespace genericga

#endif  // _GENERICGA_MULTIPOP_GA_H_
