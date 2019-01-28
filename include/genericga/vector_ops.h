#ifndef _GENERICGA_VECTOR_OPS_H_
#define _GENERICGA_VECTOR_OPS_H_

#include <algorithm>
#include <functional>
#include <numeric>
#include <type_traits>
#include <vector>

namespace genericga {
template <class T>
std::vector<int> GetOrderings(const std::vector<T>& vec);
template <class T>
std::vector<int> GetRankings(const std::vector<T>& vec);

int MinRank(int low, int high, int cur);
int MaxRank(int low, int high, int cur);
float AverageRank(int low, int high, int cur);
int CurrentRank(int low, int high, int cur);
template <class T, class Tiebreaker = std::function<int(int, int, int)>>
auto GetRankingsWithTies(const std::vector<T>& vec,
                         Tiebreaker tiebreaker = CurrentRank)
    -> std::vector<typename std::result_of<Tiebreaker(int, int, int)>::type>;
template <class T, class Set>
std::vector<T> KeyIntersection(const std::vector<T>& keys, const Set& set);
template <class T, class Set>
std::vector<T> KeyDifference(const std::vector<T>& keys, const Set& set);

// Get ordering such that ordering[n] is the index of the sorted nth element
// of vec.  In other words, vec[ordering[n]] is the nth smallest element.
// Ex. let v=[8,3,7,2,5].  v has ordering: [3,1,4,2,0].  v[3] would be index
// 0 when sorted, v[1] would be index 1 when sorted, etc.
template <class T>
std::vector<int> GetOrderings(const std::vector<T>& vec) {
  std::vector<int> indices(vec.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(),
            [&vec](int i1, int i2) { return vec[i1] < vec[i2]; });
  return indices;
}

// Get ranking such that ranking[i] is the index vec[i] would be sorted into.
// Equal values are not guaranteed a specific value.  Ex. let v=[8,3,7,2,5].  v
// has rankings [4,1,3,0,2].  v[0] would be in position 4 when sorted, v[1]
// would be in position 1, etc.
template <class T>
std::vector<int> GetRankings(const std::vector<T>& vec) {
  return GetOrderings(GetOrderings(vec));
}

// Gets ranking where ties are assigned the average of the ranks at that value.
// Ex. values of 5,6,5,1,5 have ranks: 2,4,2,0,2 (2 = mean(1,2,3)). Decimal
// ranks possible.
template <class T, class Tiebreaker>
auto GetRankingsWithTies(const std::vector<T>& vec, Tiebreaker tiebreaker)
    -> std::vector<typename std::result_of<Tiebreaker(int, int, int)>::type> {
  using Out = typename std::result_of<Tiebreaker(int, int, int)>::type;
  std::vector<Out> rankings(vec.size());
  std::vector<int> ordering = GetOrderings(vec);
  auto it = ordering.begin();
  auto end = ordering.end();
  int low_rank = 0;
  while (it != end) {
    T ref_val = vec[*it];
    auto next_it = std::find_if(
        it + 1, end, [&vec, ref_val](int i) { return vec[i] != ref_val; });
    int n_equal = std::distance(it, next_it);
    for (int cur = low_rank; it != next_it; ++it, ++cur) {
      rankings[*it] = tiebreaker(low_rank, low_rank + n_equal - 1, cur);
    }
    low_rank = low_rank + n_equal;
  }
  return rankings;
}

template <class T, class Set>
std::vector<T> KeyIntersection(const std::vector<T>& keys, const Set& set) {
  std::vector<T> int_keys;
  for (const auto& key : keys) {
    if (set.find(key) != set.end()) {
      int_keys.push_back(key);
    }
  }
  return int_keys;
}

template <class T, class Set>
std::vector<T> KeyDifference(const std::vector<T>& keys, const Set& set) {
  std::vector<T> diff_keys;
  for (const auto& key : keys) {
    if (set.find(key) == set.end()) {
      diff_keys.push_back(key);
    }
  }
  return diff_keys;
}

}  // namespace genericga

#endif  // _GENERICGA_VECTOR_OPS_H_
