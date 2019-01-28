#include "genericga/binary/byte_array_genotype.h"
#include "genericga/binary/encoding.h"

#include <boost/functional/hash.hpp>
#include <cassert>

namespace genericga {
namespace binary {

std::vector<float> ByteArrayGenotype::ToFloatArray(
    std::vector<Encoding> encodings) const {
  int bit = 0;
  std::vector<float> res;
  res.reserve(encodings.size());
  for (auto encoding : encodings) {
    unsigned int int_conversion;
    if (encoding.is_gray_coded) {
      int_conversion = FromGrayCodeBits(bit, encoding.bit_precision);
    } else {
      int_conversion = FromBits(bit, encoding.bit_precision);
    }
    bit += encoding.bit_precision;
    res.push_back((encoding.max - encoding.min) /
                      std::pow(2, encoding.bit_precision) * int_conversion +
                  encoding.min);
  }
  return res;
}

unsigned int ByteArrayGenotype::FromBits(int start_bit, int n_bits) const {
  assert(n_bits <= (sizeof(int) * CHAR_BIT));
  assert(n_bits >= 0);
  assert(start_bit >= 0);
  if (n_bits == 0) {
    return 0u;
  }
  ByteBitCoordinate start(start_bit);
  ByteBitCoordinate last(start_bit + n_bits - 1);
  assert(last.byte < data_.size());
  if (start.byte == last.byte) {
    return FromSubByte(data_[start.byte], start.bit, n_bits);
  }
  unsigned int val = 0;
  // handle last byte by shifting all unused bits off the right side
  val |= data_[last.byte] >> (CHAR_BIT - last.bit - 1);

  int n_full_bytes = (last.byte - start.byte) - 1;
  for (int i = 1; i <= n_full_bytes; ++i) {
    val |= data_[last.byte - i] << (CHAR_BIT * (i - 1) + last.bit + 1);
  }
  // handle first bit by masking off all unused bits on the left side
  val |= (GenerateMask(CHAR_BIT - start.bit) & data_[start.byte])
         << (CHAR_BIT * n_full_bytes + last.bit + 1);
  return val;
}

unsigned int ByteArrayGenotype::FromGrayCodeBits(int start_bit,
                                                 int n_bits) const {
  unsigned int val = FromBits(start_bit, n_bits);
  return GrayToBinary(val);
}

void ByteArrayGenotype::Flip(int bit_index) {
  assert(bit_index < data_.size() * CHAR_BIT);
  assert(bit_index >= 0);
  ByteBitCoordinate loc(bit_index);
  FlipBit(data_[loc.byte], loc.bit);
}

bool operator<(const ByteArrayGenotype& g1, const ByteArrayGenotype& g2) {
  return g1.data_ < g2.data_;
}

bool operator==(const ByteArrayGenotype& g1, const ByteArrayGenotype& g2) {
  return g1.data_ == g2.data_;
}

// Swaps a range of bits between two genotypes.  Starting from the left, swaps
// from first_bit to end_bit (exclusively).
void SwapBits(ByteArrayGenotype& gene1, ByteArrayGenotype& gene2, int first_bit,
              int n_bits) {
  assert(first_bit + n_bits <= gene1.NBits());
  assert(first_bit + n_bits <= gene2.NBits());
  assert(first_bit >= 0);
  if (n_bits < 0) {
    return;
  }

  ByteBitCoordinate first(first_bit);
  ByteBitCoordinate last(first_bit + n_bits - 1);

  if (first.byte == last.byte) {
    PartialByteSwap(gene1.data_[first.byte], gene2.data_[first.byte],
                    CHAR_BIT - first.bit);
    PartialByteSwap(gene1.data_[last.byte], gene2.data_[last.byte],
                    CHAR_BIT - last.bit - 1);
    return;
  }
  // Partially swap first_byte from first_sub_bit to end of the byte
  PartialByteSwap(gene1.data_[first.byte], gene2.data_[first.byte],
                  CHAR_BIT - first.bit);

  // Partially swap last_byte from start of the byte to last_sub_bit
  std::swap(gene1.data_[last.byte], gene2.data_[last.byte]);
  PartialByteSwap(gene1.data_[last.byte], gene2.data_[last.byte],
                  CHAR_BIT - last.bit - 1);

  // Fully swap all middle bytes
  std::swap_ranges(gene1.data_.begin() + first.byte + 1,
                   gene1.data_.begin() + last.byte,
                   gene2.data_.begin() + first.byte + 1);
  // for (int i = first.byte + 1; i < last.byte; ++i) {
  //   std::swap(gene1.data_[i], gene2.data_[i]);
  // }
}

// Parses all bits from start_bit to end_bit (exclusively) from a byte.
// Ex. FromSubByte(0b11001110, 2, 4) -> xx0011xx -> 0b0011 -> 3
unsigned char FromSubByte(unsigned char byte, int start_bit, int n_bits) {
  return (GenerateMask(CHAR_BIT - start_bit) & byte) >>
         (CHAR_BIT - n_bits - start_bit);
}

// From Wikipedia:
// A more efficient version for Gray codes 32 bits or fewer through the use of
// SWAR (SIMD within a register) techniques. It implements a parallel prefix
// XOR function.
unsigned int GrayToBinary(unsigned int num) {
  num = num ^ (num >> 16);
  num = num ^ (num >> 8);
  num = num ^ (num >> 4);
  num = num ^ (num >> 2);
  num = num ^ (num >> 1);
  return num;
}

unsigned char GenerateMask(const int n_bits) {
  if (n_bits == 0) return 0;
  const unsigned char src = 1u << (n_bits - 1);
  return (src - 1) ^ src;
}

void FlipBit(unsigned char& byte, int bit) {
  assert(bit >= 0);
  byte ^= 1u << (CHAR_BIT - bit - 1);
}

// XOR swap algorithm with a mask so it only swaps partially.
void PartialByteSwap(unsigned char& t1, unsigned char& t2, int n_bits) {
  assert(n_bits <= CHAR_BIT);
  if (n_bits <= 0) {
    return;
  }
  unsigned char mask = GenerateMask(n_bits);
  t1 = t1 ^ (mask & t2);
  t2 = t2 ^ (mask & t1);
  t1 = t1 ^ (mask & t2);
}

ByteBitCoordinate NextBit(ByteBitCoordinate coord) {
  return (coord.bit == CHAR_BIT - 1)
             ? ByteBitCoordinate{coord.byte + 1, 0}
             : ByteBitCoordinate{coord.byte, coord.bit + 1};
}

std::size_t hash_value(const ByteArrayGenotype& s) {
  boost::hash<std::vector<unsigned char>> hasher;
  return hasher(s.data_);
}

}  // namespace binary
}  // namespace genericga
