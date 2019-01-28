#include "genericga/vector_ops.h"

namespace genericga {

int MinRank(int low, int high, int cur) { return low; }
int MaxRank(int low, int high, int cur) { return high; }
float AverageRank(int low, int high, int cur) { return (low + high) / 2.0; }
int CurrentRank(int low, int high, int cur) { return cur; }

}  // namespace genericga
