#ifndef NUMERICALDISTS_ORDER_STATISTIC_OPS_H_
#define NUMERICALDISTS_ORDER_STATISTIC_OPS_H_

#include <eigen3/Eigen/Core>

namespace numericaldists {

Eigen::ArrayXd KthLowestOrderStatisticCDF(
    const Eigen::ArrayXd& xs, const std::vector<Eigen::ArrayXd>& cdfs, int k);
Eigen::ArrayXd KthLowestOrderStatisticPDF(
    const Eigen::ArrayXd& xs, const std::vector<Eigen::ArrayXd>& cdfs, int k);
Eigen::ArrayXd HighestOrderStatisticCDF(
    const Eigen::ArrayXd& xs, const std::vector<Eigen::ArrayXd>& cdfs);
Eigen::ArrayXd KthLowestOrderStatisticPDF(const Eigen::ArrayXd& xs,
                                          const Eigen::ArrayXd& cdf,
                                          int n_draws, int k);
Eigen::ArrayXd KthLowestOrderStatisticCDF(const Eigen::ArrayXd& xs,
                                          const Eigen::ArrayXd& cdf,
                                          int n_draws, int k);

Eigen::ArrayXd HighestOrderStatisticCDF(const Eigen::ArrayXd& xs,
                                        const Eigen::ArrayXd& cdf, int n_draws);
Eigen::ArrayXXd JointOrderStatisticPDF(const Eigen::ArrayXd& xs,
                                       const Eigen::ArrayXd& cdf, int n_draws,
                                       int j, int k);
Eigen::ArrayXXd LowestHighestJointOrderStatisticPDF(const Eigen::ArrayXd& xs,
                                                    const Eigen::ArrayXd& cdf,
                                                    int n_draws);
Eigen::ArrayXXd LowestHighestJointOrderStatisticCDF(const Eigen::ArrayXd& xs,
                                                    const Eigen::ArrayXd& cdf,
                                                    int n_draws);
}  // namespace numericaldists

#endif  // NUMERICALDISTS_ORDER_STATISTIC_OPS_H_
