#ifndef BIDDING_PLAYER_BID_H_
#define BIDDING_PLAYER_BID_H_

#include "bidding/piecewise_linear_function.h"

namespace bidding {

struct PlayerBid {
  int id;
  PiecewiseLinearFunction bid;
};

}  // namespace bidding

#endif  // BIDDING_PLAYER_BID_H_
