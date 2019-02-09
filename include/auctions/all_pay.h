#ifndef _AUCTIONS_ALL_PAY_H_
#define _AUCTIONS_ALL_PAY_H_

#include <functional>
#include <vector>

#include "numericaldists/piecewise_linear.h"

namespace auctions {

class AllPay {
 public:
  AllPay(std::vector<float> values);
  float GetValue(float bid, float value) const;
  void AcceptStrategy(std::function<float(float)> cdf, int id);
  float GetFitness(const std::function<float(float)>& cdf, int id) const;
  std::vector<float> GetFitness(
      const std::vector<numericaldists::PiecewiseLinear>& funcs, int id) const {
    std::vector<float> fits;
    fits.reserve(funcs.size());
    for (const auto& func : funcs) {
      fits.push_back(GetFitness(func, id));
    }
    return fits;
  }

 private:
  float GetIntegrand(const std::function<float(float)>& cdf, int id,
                     float bid) const;
  std::vector<std::function<float(float)>> bid_cdf_funcs_;
  std::vector<std::function<float(float)>> bid_pdf_funcs_;
  std::vector<float> values_;
  int n_players_;
};

}  // namespace auctions

#endif  // _AUCTIONS_ALL_PAY_H_
