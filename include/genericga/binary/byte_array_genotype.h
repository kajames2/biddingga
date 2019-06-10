#ifndef GENERICGA_BINARY_BYTE_ARRAY_GENOTYPE_H_
#define GENERICGA_BINARY_BYTE_ARRAY_GENOTYPE_H_

#include <eigen3/Eigen/Core>

#include <cassert>
#include <climits>
#include <vector>

#include "genericga/binary/float_encoding.h"
#include "genericga/binary/int_encoding.h"

namespace genericga {
namespace binary {

struct ByteBitCoordinate;

// Genotype representation as a byte array that can be flipped and swapped with
// other genotypes.  Stores everything in little-endian at both byte and bit
// level (bit 0 is left-most bit), and bit N is rightmost bit).  The endianness
// should be considered an implementation detail.
class ByteArrayGenotype {
 public:
  ByteArrayGenotype() {}
  explicit ByteArrayGenotype(std::vector<unsigned char> data) : data_(data) {}
  void Flip(int bit);
  std::vector<int> ToIntArray(std::vector<IntEncoding> encodings) const;
  std::vector<float> ToFloatArray(std::vector<FloatEncoding> encodings) const;
  Eigen::ArrayXf ToEigenFloatArray(std::vector<FloatEncoding> encodings) const;
  unsigned int FromGrayCodeBits(int start_bit, int n_bits) const;
  unsigned int FromBits(int start_bit, int n_bits) const;
  friend void SwapBits(ByteArrayGenotype& gene1, ByteArrayGenotype& gene2,
                       int start_bit, int n_bits);
  unsigned char& operator[](int index) { return data_[index]; }
  int NBytes() const { return data_.size(); }
  int NBits() const { return data_.size() * CHAR_BIT; }
  friend bool operator<(const ByteArrayGenotype& g1,
                        const ByteArrayGenotype& g2);
  friend bool operator==(const ByteArrayGenotype& g1,
                         const ByteArrayGenotype& g2);
  friend std::size_t hash_value(const ByteArrayGenotype& s);

 private:
  std::vector<unsigned char> data_;
  unsigned char GetBit(ByteBitCoordinate coord);
};

unsigned char FromSubByte(unsigned char byte, int start_bit, int n_bits);
unsigned int GrayToBinary(unsigned int num);
// Generates a 1-byte mask.  n_bits is the number of 1's (from the right).
// Ex. GenerateMask(5) -> 00011111.
unsigned char GenerateMask(const int n_bits);
void FlipBit(unsigned char& byte, int bit);
void PartialByteSwap(unsigned char& t1, unsigned char& t2, int n_bits);
void SwapBits(ByteArrayGenotype& gene1, ByteArrayGenotype& gene2, int first_bit,
              int n_bits);

}  // namespace binary
}  // namespace genericga

#endif  // GENERICGA_BINARY_BYTE_ARRAY_GENOTYPE_H_
