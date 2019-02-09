#ifndef _GENERICGA_SELECTOR_ELITISM_DECORATOR_H_
#define _GENERICGA_SELECTOR_ELITISM_DECORATOR_H_

#include <memory>
#include <vector>

#include "genericga/selector.h"
#include "genericga/selector/keep_best.h"

namespace genericga {
namespace selector {

class ElitismDecorator : public Selector {
 public:
  ElitismDecorator(std::unique_ptr<Selector> sel, int n_elites);
  std::vector<int> SelectIndices(const std::vector<float>& fitnesses,
                                 const std::vector<int>& counts,
                                 int n) override;

 private:
  std::unique_ptr<Selector> sel_;
  std::unique_ptr<KeepBest> elite_sel_;
  int n_elites_;
};

}  // namespace selector
}  // namespace genericga

#endif  // _GENERICGA_SELECTOR_ELITISM_DECORATOR_H_
