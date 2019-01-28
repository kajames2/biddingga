#ifndef _SAMPLE_PHENOTYPE_CONVERTER_H_
#define _SAMPLE_PHENOTYPE_CONVERTER_H_

namespace gatests {

class SamplePhenotypeConverter {
public:
  explicit SamplePhenotypeConverter(int n) : n_(n) {}
  int operator()(const int &genotype) const { return n_ - genotype * genotype; }

private:
  int n_;
};

} // namespace gatests

#endif // _SAMPLE_PHENOTYPE_CONVERTER_H_
