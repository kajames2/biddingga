#include "numericaldists/combination_generation.h"

#include <algorithm>
#include <functional>
#include <vector>

namespace numericaldists {

unsigned int GetFirstCanonicalCombination(int r) {
  if (r == 0) return 0;
  const unsigned char src = 1u << (r - 1);
  return (src - 1) ^ src;
}

// Algorithm found at: https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
unsigned int GetNextCanonicalCombination(unsigned int v) {
  unsigned int t = v | (v - 1); // t gets v's least significant 0 bits set to 1
  // Next set to 1 the most significant bit to change,
  // set to 0 the least significant ones, and add the necessary 1 bits.
  return (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(v) + 1));
}

}  // namespace numericaldists
