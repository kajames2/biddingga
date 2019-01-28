#include <gtest/gtest.h>

#include <vector>

#include "genericga/vector_ops.h"

namespace gatests {

class VectorOpsTest : public ::testing::Test {
public:
  VectorOpsTest() {}

protected:
  virtual void SetUp() {
    vec.push_back(4);
    vec.push_back(0);
    vec.push_back(-2);
    vec.push_back(5);
    vec.push_back(4);
    vec.push_back(1);

  }
  std::vector<int> vec;
};

TEST_F(VectorOpsTest, GetOrderingsTest) {
  auto rankings = genericga::GetOrderings(vec);
  ASSERT_EQ(1, rankings[1]);
  ASSERT_EQ(5, rankings[2]);
}

TEST_F(VectorOpsTest, GetRankingsTest) {
  auto rankings = genericga::GetRankings(vec);
  ASSERT_EQ(1, rankings[1]);
  ASSERT_EQ(0, rankings[2]);
}

TEST_F(VectorOpsTest, GetRankingsWithTiesTest) {
  auto rankings = genericga::GetRankingsWithTies(vec);
  ASSERT_EQ(0, rankings[2]);
  ASSERT_EQ(3, rankings[0]);
}

TEST_F(VectorOpsTest, GetRankingsWithTiesFloatTest) {
  auto rankings = genericga::GetRankingsWithTies(vec, genericga::AverageRank);
  ASSERT_FLOAT_EQ(0, rankings[2]);
  ASSERT_FLOAT_EQ(3.5, rankings[0]);
}

TEST_F(VectorOpsTest, KeyDifferenceTest) {
  std::vector<int> vec1 = {1,2,9,4};
  std::set<int> ref = {1,8,9};
  std::vector<int> diff = {2,4};
  EXPECT_EQ(diff, genericga::KeyDifference(vec1, ref));
}

TEST_F(VectorOpsTest, KeyIntersectionTest) {
  std::vector<int> vec1 = {1,2,9,4};
  std::set<int> ref = {1,8,9};
  std::vector<int> in = {1,9};
  EXPECT_EQ(in, genericga::KeyIntersection(vec1, ref));
}

} // namespace gatests
