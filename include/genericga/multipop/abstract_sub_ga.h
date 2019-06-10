#ifndef GENERICGA_MULTIPOP_ABSTRACT_SUB_GA_H_
#define GENERICGA_MULTIPOP_ABSTRACT_SUB_GA_H_

#include <vector>

namespace genericga {
namespace multipop {

template <class Environment>
class AbstractSubGA {
 public:
  virtual void RunRound(const std::vector<Environment*>& envs) = 0;
  virtual void SubmitPlayStrat(Environment& env) = 0;
  virtual ~AbstractSubGA() {}
  
  void SubmitPlayStrats(std::vector<Environment*>& envs) {
    for (auto env : envs) {
      SubmitPlayStrat(env);
    }
  }
  int GetID() const { return id_; }
  int GetPriority() const { return priority_; }
  void SetPriority(int priority) { priority_ = priority; }

 protected:
    explicit AbstractSubGA(int id, int priority = 0)
      : id_(id), priority_(priority) {}
 private:
  int id_;
  int priority_;
};

}  // namespace multipop
}  // namespace genericga

#endif  // GENERICGA_MULTIPOP_ABSTRACT_SUB_GA_H_
