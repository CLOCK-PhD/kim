#include "file_reader.h"

#include <fstream>
#include <iostream>

using namespace std;

FileReader::FileReader(const char *filename, size_t k):
  filename(filename), ifs(filename),
  line(0), col(0),
  k(0), read_id(), kmer(), kmer_pos(-1) {
}  

FileReader::~FileReader() {
  ifs.close();
}

const string &FileReader::getNextKmer() {
  // TODO
  cerr << "TODO:" << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << endl;
  return kmer;
}

const string &FileReader::getCurrentKmer() const {
  return kmer;
}

const string &FileReader::getCurrentRead() const {
  return read_id;
}

size_t FileReader::getCurrentKmerRelativePosition() const {
  return kmer_pos;
}
