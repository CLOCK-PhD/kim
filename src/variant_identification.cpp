#include "variant_identification.h"

#include <algorithm>
#include <iostream>

using namespace std;

typedef unordered_multimap<VariantIdentification::ReadID_type, size_t> _mymap;

VariantIdentification::VariantIdentification(): informations() {
}

VariantIdentification &VariantIdentification::add(VariantIdentification::ReadID_type read_id, size_t pos) {
  informations.emplace(read_id, pos);
  return *this;
}

bool cmp_key(_mymap::value_type i, _mymap::value_type j) {
  return i.first == j.first;
}

list<VariantIdentification::ReadID_type> VariantIdentification::getReads() const {
  list<VariantIdentification::ReadID_type> l;
  _mymap::const_iterator it_prev = informations.cbegin();
  for (_mymap::const_iterator it = informations.cbegin(); it != informations.cend(); ++it) {
    if (it_prev->first != it->first) {
      l.push_back(it->first);
      it_prev = it;
    }
  }
  if (it_prev != informations.cend()) {
    l.push_back(it_prev->first);
  }
  return l;
}

double VariantIdentification::getReadScore(VariantIdentification::ReadID_type read_id) const {
  double v = 0;
  pair<_mymap::const_iterator, _mymap::const_iterator> range = informations.equal_range(read_id);
  for (_mymap::const_iterator it = range.first; it != range.second; ++it) {
    // TODO with it->second (one of the position)
    ++v;
  }
  cerr << "TODO:" << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << endl;
  return v;
}
