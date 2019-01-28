#ifndef _BIDDING_PLAYER_BID_H_
#define _BIDDING_PLAYER_BID_H_

#include "bidding/piecewise_linear_function.h"

namespace bidding {

struct PlayerBid {
  int id;
  PiecewiseLinearFunction bid;
};

}  // namespace bidding

#endif  // _BIDDING_PLAYER_BID_H_
