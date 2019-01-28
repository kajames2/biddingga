#ifndef _GENERICGA_BINARY_GENOTYPE_H_
#define _GENERICGA_BINARY_GENOTYPE_H_

namespace genericga {
namespace binary {

class Genotype {
 public:
  virtual void Flip(int bit) = 0;
  virtual unsigned int FromBits(int start_bit, int n_bits) const = 0;
  unsigned int FromGrayCodeBits(int start_bit, int n_bits) const;
  virtual void SwapBits(Genotype& gene2, int start_bit, int n_bits);
  virtual int NBits() const = 0;
  virtual ~Genotype() {}
};

}  // namespace binary
}  // namespace genericga

#endif  // _GENERICGA_BINARY_GENOTYPE_H_
