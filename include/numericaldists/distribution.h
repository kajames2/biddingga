#ifndef _NUMERICALDISTS_DISTRIBUTION_H_
#define _NUMERICALDISTS_DISTRIBUTION_H_

#include <boost/math/distributions.hpp>

#include <memory>

namespace numericaldists {

// Code adapted from Sean Parent's talk on runtime polymorphism.
// Is an example of type erasure intended to wrap boost's distributions so that
// I can put different distributions into a vector and use any combination of
// value distributions in auctions.
class Distribution {
 public:
  template <typename T>
  Distribution(T x) : self_(std::make_shared<model<T>>(std::move(x))) {}

  friend float cdf(const Distribution& dist, float x) {
    return dist.self_->cdf_(x);
  }
  friend float pdf(const Distribution& dist, float x) {
    return dist.self_->pdf_(x);
  }
  friend float quantile(const Distribution& dist, float x) {
    return dist.self_->quantile_(x);
  }

 private:
  struct DistributionConcept {
    virtual ~DistributionConcept() = default;
    virtual float cdf_(float x) const = 0;
    virtual float pdf_(float x) const = 0;
    virtual float quantile_(float x) const = 0;
  };
  template <typename T>
  struct model final : DistributionConcept {
    model(T x) : dist_(std::move(x)) {}
    float cdf_(float x) const override { return boost::math::cdf(dist_, x); }
    float pdf_(float x) const override { return boost::math::pdf(dist_, x); }
    float quantile_(float x) const override {
      return boost::math::quantile(dist_, x);
    }
    T dist_;
  };

  std::shared_ptr<const DistributionConcept> self_;
};

float lower(const Distribution& dist);
float upper(const Distribution& dist);

float lower(const std::vector<Distribution>& dists);
float upper(const std::vector<Distribution>& dists);

}  // namespace numericaldists

#endif  // _NUMERICALDISTS_DISTRIBUTION_H_
