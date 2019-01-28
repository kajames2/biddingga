#include <gtest/gtest.h>

#include <vector>

#include "numericaldists/combination_generation.h"

namespace gatests {

using namespace numericaldists;

class CombinationGenerationTest : public ::testing::Test {
 public:
  CombinationGenerationTest() {}

 protected:
  virtual void SetUp() {}
};

TEST_F(CombinationGenerationTest, GenerateFirstsTest) {
  EXPECT_EQ(0, GetFirstCanonicalCombination(0));
  EXPECT_EQ(7, GetFirstCanonicalCombination(3));
  EXPECT_EQ(31, GetFirstCanonicalCombination(5));
}

TEST_F(CombinationGenerationTest, GetNextTest) {
  unsigned int val = GetFirstCanonicalCombination(2);
  EXPECT_EQ(3, val);
  val = GetNextCanonicalCombination(val);
  EXPECT_EQ(5, val);
  val = GetNextCanonicalCombination(val);
  EXPECT_EQ(6, val);
}
}  // namespace gatests
