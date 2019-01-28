#ifndef _GENERICGA_MULTIPOP_ABSTRACT_SUB_GA_H_
#define _GENERICGA_MULTIPOP_ABSTRACT_SUB_GA_H_

#include <vector>

namespace genericga {
namespace multipop {

template <class Environment>
class AbstractSubGA {
 public:
  AbstractSubGA(int id, int priority = 0) : id_(id), priority_(priority) {}
  void SubmitPlayStrats(std::vector<Environment>& envs) {
    for (auto& env : envs) {
      SubmitPlayStrat(env);
    }
  }
  virtual void RunRound(const std::vector<Environment>& envs);
  int GetID() const { return id_; }
  int GetPriority() const { return priority_; }
  void SetPriority(int priority) { priority_ = priority; }
  virtual void SubmitPlayStrat(Environment& env);
  virtual ~AbstractSubGA() {}
  
 private:
  int id_;
  int priority_;
};

}  // namespace multipop
}  // namespace genericga

#endif  // _GENERICGA_MULTIPOP_ABSTRACT_SUB_GA_H_
