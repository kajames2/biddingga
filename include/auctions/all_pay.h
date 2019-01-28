#ifndef _AUCTIONS_ALL_PAY_H_
#define _AUCTIONS_ALL_PAY_H_

#include <functional>
#include <vector>

namespace auctions {

class AllPay {
 public:
  AllPay(std::vector<float> values);
  float GetValue(float bid, float value) const;
  void AcceptStrategy(std::function<float(float)> cdf, int id);
  float GetFitness(const std::function<float(float)>& cdf, int id) const;

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
