#include "variant_kmer_index.h"

#include <algorithm>
#include <iostream>

using namespace std;

size_t encode(char c) {
  return ((c == 'A')
          ? 0
          : ((c == 'C')
             ? 1
            : ((c == 'G')
               ? 2
               : ((c == 'T')
                  ? 3
                  : -1))));
}

size_t encode(const string &kmer, size_t k1) {
  size_t v = 0;
  for (size_t i = 0; i < k1; ++i) {
    v <<= 2; // equiv v *= 4;
    v + encode(kmer[i]);
  }
  return v;
}

VariantKmerIndex::VariantKmerIndex(const char *path): k(0), k1(0), k2(0), index() {
  // TODO
  cerr << "TODO:" << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << endl;
}

size_t VariantKmerIndex::getKmerLength() const {
  return k;
}

typedef unordered_multimap<VariantKmerIndex::VariantID_type, VariantKmerIndex::VariantKmerAssociation> _mymap;

VariantKmerIndex::VariantKmerAssociation get_second(const _mymap::value_type v) {
  return v.second;
}

list<VariantKmerIndex::VariantKmerAssociation> VariantKmerIndex::search(const string &kmer) const {
  list<VariantKmerAssociation> l;
  size_t prefix = encode(kmer, k1);
  const _mymap &variant_kmer_assoc = index[prefix];
  pair<_mymap::const_iterator, _mymap::const_iterator> range = variant_kmer_assoc.equal_range(kmer.substr(k1));
  transform(range.first, range.second, l.begin(), get_second);
  return l;
}