#ifndef _GENERICGA_CHILDREN_FACTORY_H_
#define _GENERICGA_CHILDREN_FACTORY_H_

#include <memory>
#include <vector>
#include <iostream>
#include "genericga/crossover.h"
#include "genericga/genotype_population.h"
#include "genericga/mutator.h"
#include "genericga/selector.h"
#include "genericga/selector/tournament.h"

namespace genericga {

// Generates new genotypes (children) by performing selection, crossover, and
// mutation on a population.
template <class Gen>
class ChildrenFactory {
 public:
  ChildrenFactory(std::unique_ptr<Crossover<Gen>> crossover,
                  std::unique_ptr<Mutator<Gen>> mutator,
                  std::unique_ptr<Selector> parent_selector = std::make_unique<selector::Tournament>(2));
  std::vector<Gen> GetChildren(const GenotypePopulation<Gen>& pop_,
                               int n_children_);

 private:
  // Pairs up adjacent parents, and conducts crossover.  If odd # of parents,
  // last parent is untouched.
  void ConductCrossover(std::vector<Gen>& children);
  void ConductMutation(std::vector<Gen>& children);

  std::unique_ptr<Crossover<Gen>> crossover_;
  std::unique_ptr<Mutator<Gen>> mutator_;
  std::unique_ptr<Selector> parent_selector_;
};

template <class Gen>
ChildrenFactory<Gen>::ChildrenFactory(std::unique_ptr<Crossover<Gen>> crossover,
                                      std::unique_ptr<Mutator<Gen>> mutator,
                                      std::unique_ptr<Selector> parent_selector)
    : crossover_(std::move(crossover)),
      mutator_(std::move(mutator)),
      parent_selector_(std::move(parent_selector)) {}

template <class Gen>
std::vector<Gen> ChildrenFactory<Gen>::GetChildren(
    const GenotypePopulation<Gen>& pop, int n_children) {
  auto children = pop.SelectGenotypes(*parent_selector_, n_children);
  ConductCrossover(children);
  ConductMutation(children);
  return children;
}

template <class Gen>
void ChildrenFactory<Gen>::ConductCrossover(std::vector<Gen>& children) {
  // Gotta stop one before the end since we use i+1 as the neighboring
  // geneotype.
  for (int i = 0; i < children.size() - 1; i += 2) {
    (*crossover_)(children[i], children[i + 1]);
  }
}

template <class Gen>
void ChildrenFactory<Gen>::ConductMutation(std::vector<Gen>& children) {
  for (auto& child : children) {
    (*mutator_)(child);
  }
}

}  // namespace genericga

#endif  // _GENERICGA_CHILDREN_FACTORY_H_
