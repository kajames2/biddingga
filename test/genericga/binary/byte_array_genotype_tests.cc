#include <gtest/gtest.h>

#include <vector>

#include "genericga/binary/byte_array_genotype.h"

namespace gatests {

class ByteArrayGenotypeTest : public ::testing::Test {
 public:
  ByteArrayGenotypeTest() {}

 protected:
  virtual void SetUp() {}

  genericga::binary::ByteArrayGenotype gene1 = genericga::binary::ByteArrayGenotype(
      std::vector<unsigned char>{0xFF, 0x00, 0xFF, 0x00});
  genericga::binary::ByteArrayGenotype gene2 = genericga::binary::ByteArrayGenotype(
      std::vector<unsigned char>{0x00, 0xFF, 0x00, 0xFF});
};

typedef ByteArrayGenotypeTest ByteArrayGenotypeDeathTest;
using byte = unsigned char;

TEST_F(ByteArrayGenotypeTest, PartialByteSwapTest) {
  byte a = static_cast<byte>(0x00);
  byte b = static_cast<byte>(0xFF);
  genericga::binary::PartialByteSwap(a, b, 3);
  EXPECT_EQ(a, static_cast<byte>(0x07));
  EXPECT_EQ(b, static_cast<byte>(0xF8));
}

TEST_F(ByteArrayGenotypeTest, PartialByteSwapNoTest) {
  byte a = static_cast<byte>(0x00);
  byte b = static_cast<byte>(0xFF);
  genericga::binary::PartialByteSwap(a, b, 0);
  EXPECT_EQ(a, static_cast<byte>(0x00));
  EXPECT_EQ(b, static_cast<byte>(0xFF));
}

TEST_F(ByteArrayGenotypeTest, PartialByteSwapFullTest) {
  byte a = static_cast<byte>(0x00);
  byte b = static_cast<byte>(0xFF);
  genericga::binary::PartialByteSwap(a, b, 8);
  EXPECT_EQ(a, static_cast<byte>(0xFF));
  EXPECT_EQ(b, static_cast<byte>(0x00));
}

TEST_F(ByteArrayGenotypeDeathTest, PartialByteSwapFullDeathTest) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  byte a = static_cast<byte>(0x00);
  byte b = static_cast<byte>(0xFF);
  EXPECT_DEATH(genericga::binary::PartialByteSwap(a, b, 9), "");
}

TEST_F(ByteArrayGenotypeTest, GrayToBinaryTest) {
  byte a = 0x08;
  EXPECT_EQ(15, genericga::binary::GrayToBinary(a));
}

TEST_F(ByteArrayGenotypeTest, FromBitsTest) {
  EXPECT_EQ(0u, gene1.FromBits(4, 0));
  EXPECT_EQ(3u, gene1.FromBits(6, 2));
  EXPECT_EQ(static_cast<unsigned int>(0x601), gene1.FromBits(6, 11));
}

TEST_F(ByteArrayGenotypeTest, FlipTest) {
  gene1.Flip(0);
  EXPECT_EQ(127u, gene1.FromBits(0, 8));
}

TEST_F(ByteArrayGenotypeTest, SwapBitsFullBytesTest) {
  genericga::binary::SwapBits(gene1, gene2, 6, 12);
  EXPECT_EQ(static_cast<unsigned int>(0xFCFF3F00), gene1.FromBits(0, 32));
}

TEST_F(ByteArrayGenotypeTest, SwapBitsOffsetByteTest) {
  genericga::binary::SwapBits(gene1, gene2, 6, 8);
  EXPECT_EQ(static_cast<unsigned int>(0xFCFCFF00), gene1.FromBits(0, 32));
}

TEST_F(ByteArrayGenotypeTest, SwapBitsSubByteTest) {
  genericga::binary::SwapBits(gene1, gene2, 2, 3);
  EXPECT_EQ(static_cast<unsigned int>(0xC700FF00), gene1.FromBits(0, 32));
}

TEST_F(ByteArrayGenotypeTest, SwapBitsFullTest) {
  genericga::binary::SwapBits(gene1, gene2, 0, 32);
  EXPECT_EQ(static_cast<unsigned int>(0x00FF00FF), gene1.FromBits(0, 32));
}

TEST_F(ByteArrayGenotypeTest, SwapBitsNoneTest) {
  genericga::binary::SwapBits(gene1, gene2, 0, 0);
  EXPECT_EQ(static_cast<unsigned int>(0xFF00FF00), gene1.FromBits(0, 32));
}
}  // namespace gatests
