/******************************************************************************
*                                                                             *
*  Copyright © 2024      -- IGH / LIRMM / CNRS / UM                           *
*                           (Institut de Génétique Humaine /                  *
*                           Laboratoire d'Informatique, de Robotique et de    *
*                           Microélectronique de Montpellier /                *
*                           Centre National de la Recherche Scientifique /    *
*                           Université de Montpellier)                        *
*                                                                             *
*                                                                             *
*  Auteurs/Authors:                                                           *
*    - Rémy COSTA       <remy.costa@igh.cnrs.fr>                              *
*    - William RITCHIE  <william.ritchie@igh.cnrs.fr>                         *
*    - Alban MANCHERON  <alban.mancheron@lirmm.fr>                            *
*                                                                             *
*                                                                             *
*  Programmeurs/Programmers:                                                  *
*    - Rémy COSTA       <remy.costa@igh.cnrs.fr>                              *
*    - Alban MANCHERON  <alban.mancheron@lirmm.fr>                            *
*                                                                             *
*                                                                             *
*  Contact:                                                                   *
*    - KIM list         <kim@lirmm.fr>                                        *
*                                                                             *
*  -------------------------------------------------------------------------  *
*                                                                             *
*  Ce logiciel  est un  programme informatique  permettant  d'identifier des  *
*  variations génomiques à partir de données brutes de séquençage.            *
*                                                                             *
*  Ce logiciel est régi par la  licence CeCILL  soumise au droit français et  *
*  respectant les principes  de diffusion des logiciels libres.  Vous pouvez  *
*  utiliser, modifier et/ou redistribuer ce programme sous les conditions de  *
*  la licence CeCILL telle que diffusée par  le CEA,  le CNRS et l'INRIA sur  *
*  le site "http://www.cecill.info".                                          *
*                                                                             *
*  En contrepartie de l'accessibilité au code source et des droits de copie,  *
*  de modification et de redistribution accordés par cette licence, il n'est  *
*  offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,  *
*  seule une responsabilité  restreinte pèse  sur l'auteur du programme,  le  *
*  titulaire des droits patrimoniaux et les concédants successifs.            *
*                                                                             *
*  À  cet égard  l'attention de  l'utilisateur est  attirée sur  les risques  *
*  associés au   chargement, à   l'utilisation, à  la modification  et/ou au  *
*  développement et   à la reproduction du  logiciel par l'utilisateur étant  *
*  donné  sa spécificité  de logiciel libre,  qui peut le rendre  complexe à  *
*  manipuler et qui le réserve donc à des développeurs et des professionnels  *
*  avertis  possédant  des  connaissances  informatiques  approfondies.  Les  *
*  utilisateurs sont   donc invités   à charger   et tester  l'adéquation du  *
*  logiciel à   leurs besoins  dans des  conditions permettant  d'assurer la  *
*  sécurité de leurs systèmes et ou de leurs données et, plus  généralement,  *
*  à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.         *
*                                                                             *
*  Le fait que   vous puissiez accéder à cet  en-tête signifie que vous avez  *
*  pris connaissance de la licence CeCILL,   et que vous en avez accepté les  *
*  termes.                                                                    *
*                                                                             *
*  -------------------------------------------------------------------------  *
*                                                                             *
*  This software is a computer program whose purpose is to indentify genomic  *
*  from raw sequencing data.                                                  *
*                                                                             *
*  This software is governed by the CeCILL license under French law and       *
*  abiding by the rules of distribution of free software. You can use,        *
*  modify and/ or redistribute the software under the terms of the CeCILL     *
*  license as circulated by CEA, CNRS and INRIA at the following URL          *
*  "http://www.cecill.info".                                                  *
*                                                                             *
*  As a counterpart to the access to the source code and rights to copy,      *
*  modify and redistribute granted by the license, users are provided only    *
*  with a limited warranty and the software's author, the holder of the       *
*  economic rights, and the successive licensors have only limited            *
*  liability.                                                                 *
*                                                                             *
*  In this respect, the user's attention is drawn to the risks associated     *
*  with loading, using, modifying and/or developing or reproducing the        *
*  software by the user in light of its specific status of free software,     *
*  that may mean that it is complicated to manipulate, and that also          *
*  therefore means that it is reserved for developers and experienced         *
*  professionals having in-depth computer knowledge. Users are therefore      *
*  encouraged to load and test the software's suitability as regards their    *
*  requirements in conditions enabling the security of their systems and/or   *
*  data to be ensured and, more generally, to use and operate it in the same  *
*  conditions as regards security.                                            *
*                                                                             *
*  The fact that you are presently reading this means that you have had       *
*  knowledge of the CeCILL license and that you accept its terms.             *
*                                                                             *
******************************************************************************/

#include <kim_settings.h>
#include <dna_file_index.h>
#include <dna_file_reader.h>

#include <cstdlib>
#include <iostream>
#include <vector>
#include <list>
#include <set>
#include <string>
#include <regex>
#ifdef NDEBUG
#  undef NDEBUG
#endif
#include <cassert>

using namespace std;
using namespace kim;

bool operator==(const FileReader::FileState &s1, const FileReader::FileState &s2) {
  return ((s1.filename == s2.filename)
          && (s1.line == s2.line)
          && (s1.column == s2.column)
          && (s1.pos == s2.pos));
}

bool operator!=(const FileReader::FileState &s1, const FileReader::FileState &s2) {
  return !operator==(s1, s2);
}

bool operator==(const DNAFileReader::FileState &s1, const DNAFileReader::FileState &s2) {
  return (((FileReader::FileState)s1  == (FileReader::FileState)s2)
          && (s1.format == s2.format)
          && (s1.k == s2.k)
          && (s1.start_symbol_expected == s2.start_symbol_expected)
          && (s1.nb_nucleotides == s2.nb_nucleotides)
          && (s1.current_sequence_id == s2.current_sequence_id)
          && (s1.current_sequence_description == s2.current_sequence_description)
          && (s1.current_sequence_length == s2.current_sequence_length)
          && (s1.kmer_aux == s2.kmer_aux)
          && (s1.kmer_start_pos == s2.kmer_start_pos)
          && (s1.kmer == s2.kmer)
          && (s1.sequence_start_file_state == s2.sequence_start_file_state));
}

bool operator!=(const DNAFileReader::FileState &s1, const DNAFileReader::FileState &s2) {
  return !operator==(s1, s2);
}

ostream &operator<<(ostream &os, const FileReader::FileState &s) {
  os << "'" << s.filename << ":" << s.line << ":" << s.column;
  return os;
}

ostream &operator<<(ostream &os, const DNAFileReader::FileState &s) {
  os << "DNAFile[" << s.format << "]:" << (FileReader::FileState) s << "\n"
     << "- k: " << s.k
     << ", start_symbol_expected: " << s.start_symbol_expected
     << ", nb_nucleotides: " << s.nb_nucleotides
     << ", current_sequence_id: " << s.current_sequence_id
     << ", current_sequence_description: " << s.current_sequence_description
     << ", current_sequence_length: " << s.current_sequence_length
     << ", kmer_aux: " << s.kmer_aux
     << ", kmer_start_pos: " << s.kmer_start_pos
     << ", kmer: " << s.kmer
     << ", sequence_start_file_state: " << s.sequence_start_file_state
     << "\n";
  return os;
}

ostream &operator<<(ostream &os, const DNAFileIndex::Status &status) {
  switch (status) {
  case DNAFileIndex::SUCCESS:
    os << "SUCCESS"; break;
  case DNAFileIndex::SEQUENCE_NOT_FOUND:
    os << "SEQUENCE_NOT_FOUND"; break;
  case DNAFileIndex::POSITION_NOT_FOUND:
    os << "POSITION_NOT_FOUND"; break;
  case DNAFileIndex::FILE_NOT_FOUND:
    os << "FILE_NOT_FOUND"; break;
  case DNAFileIndex::FILE_PARSE_ERROR:
    os << "FILE_PARSE_ERROR"; break;
  case DNAFileIndex::SEQUENCE_NAME_DUPLICATED:
    os << "SEQUENCE_NAME_DUPLICATED"; break;
  }
  return os;
}

class DNAFileIndex_Tester {

private:

  typedef pair<string, DNAFileReader::Format> File_t;

  DNAFileIndex _index;
  DNAFileReader _reader;

  size_t _expected_nb_sequences;
  size_t _expected_nb_nucleotides;
  list<File_t> _expected_files;

public:

  DNAFileIndex_Tester(size_t k = 5, size_t bookmark_distance = 10):
    _index(bookmark_distance), _reader(k),
    _expected_nb_sequences(0), _expected_nb_nucleotides(0),
    _expected_files() {
    cout << "* Testing DNAFileIndex with:" << endl
         << "  - k = " << _reader.k() << " (expecting " << k << ")" << endl
         << "  - bookmark distance = " << _index.bookmark_distance
         << " (expecting " << bookmark_distance << ")" << endl
         << endl;
    assert(_reader.k() == k);
    assert(_index.bookmark_distance == bookmark_distance);
    // The testing files are inconsistent and thus require consistency
    // checking.
    consistencyChecking(true);
    warn(true);
  }

  void consistencyChecking(bool check_consistency) {
    _index.consistencyChecking(check_consistency);
    _reader.check_consistency = check_consistency;
  }

  void warn(bool value) {
    _index.warn(value);
    _reader.warn = value;
  }

  void clear() {
    cout << "* Clearing index" << endl;
    _reader.close();
    _index.clear();
    _expected_nb_sequences = 0;
    _expected_nb_nucleotides = 0;
    _expected_files.clear();
  }

  string findFile(const string &fname) {
    cout << "Looking for file '" << fname << "' in the package directories." << endl;
    static const char *dirs[] = {
      "./",
      PACKAGE_DATADIR "/",
      /* The following directories are mostly for development purpose */
      "resources/",
      "../resources/",
      "../",
      SRCDIR "/../resources/",
      /* The last array entry must be NULL */
      NULL
    };
    string new_fname = FileReader::findFile(fname, dirs);
    cout << "File path should not be empty: '" << new_fname << "'" << endl;
    assert(!new_fname.empty());
    return new_fname;
  }

  void setFilename(const string &fname, DNAFileReader::Format fmt) {
    cout << "* Trying to set filename '" << fname << "'" << endl;
    string new_fname = findFile(fname);
    _reader.open(new_fname);
    assert(_reader);
    cout << "File '" << new_fname << "' is supposed to be '" << fmt << "' formatted" << endl;
    list<string> indexed_files = _index.getFilenames();
    if (find(indexed_files.begin(), indexed_files.end(), new_fname) == indexed_files.end()) {
      cout << "File '" << new_fname << "' is added to expected files" << endl;
      _expected_files.emplace_back(new_fname, fmt);
    } else {
      cout << "File '" << new_fname << "' already in index" << endl;
    }
    DNAFileReader::Format fmt2 = _reader.getFormat();
    cout << "Detected format is '" << fmt2 << "'" << endl
         << endl;
    assert(fmt == fmt2);
    _index.sync(_reader);
  }

  void check() {

    cout << "* Checking index" << endl
         << "- Index statistics:" << endl;

    size_t nb_nucleotides = 0;
    size_t nb_sequences = 0;

    size_t nb_files = _index.getNumberOfFiles();
    cout << "  - Files [" << nb_files
         << " (expecting " << _expected_files.size() << ")]:"
         << endl;
    assert(nb_files == _expected_files.size());

    list<File_t>::const_iterator it = _expected_files.begin();
    for (const auto &f: _index.getFilenames()) {
      size_t nb_nucleotides_in_file = 0;

      DNAFileReader::Format fmt = _index.getFileFormat(f);
      cout << "    - Name: '" << f << "' (expecting '" << it->first << "')" << endl
           << "      Format: '" << fmt << "' (expecting '" << it->second << "')" << endl;
      assert(f == it->first);
      assert(fmt == it->second);
      ++it;
      nb_sequences += _index.getNumberOfSequences(f);

      cout << "      Sequences [" << nb_sequences
           << " (expecting " << _index.getNumberOfSequences(f) << ")]:"
           << endl;
      for (const auto &s: _index.getSequenceNames(f)) {
        cout << "      - Name: '" << s << "'" << endl;
        size_t nb_nucleotides_in_sequence = _index.getNumberOfNucleotides(f, s);
        cout << "        Nb of total nucleotides: " << nb_nucleotides_in_sequence << endl;
        nb_nucleotides_in_file += nb_nucleotides_in_sequence;
      }

      cout << "      Nb total Nucleotides for file: " << nb_nucleotides_in_file << endl;
      nb_nucleotides += nb_nucleotides_in_file;

    }

    cout << "    Nb total Nucleotides for index: " << nb_nucleotides
         << " (expecting " << _expected_nb_nucleotides << ")" << endl;
    assert(nb_nucleotides == _expected_nb_nucleotides);

    cout << "    Number of sequences (computed): " << nb_sequences
         << " (expecting " << _expected_nb_sequences << "):"
         << endl;
    assert(nb_sequences == _expected_nb_sequences);

    nb_sequences = _index.getNumberOfSequences();
    cout << "    Number of sequences (index query): " << nb_sequences
         << " (expecting " << _expected_nb_sequences << "):"
         << endl;

    for (const auto &s: _index.getSequenceNames()) {
      const string file = _index.getSequenceFile(s);
      const size_t id = _index.getSequenceID(file, s);
      cout << "    - '" << s << "' [defined in '" << file << "' with id " << id << "]" << endl;
      const size_t id_other = _index.getSequenceID(s);
      cout << "      Retrieving ID directly by its name returns " << id_other << endl;
      assert(id == id_other);
      const string name = _index.getSequenceName(file, id);
      cout << "      Retrieving name from ID " << id << " returns '" << name << "'" << endl;
      assert(s == name);
      const string description1 = _index.getSequenceDescription(file, s);
      const string description2 = _index.getSequenceDescription(file, id);
      const string description3 = _index.getSequenceDescription(s);
      cout << "      Retrieving description from:" << endl;
      cout << "      - name in file: '" << description1 << "'" << endl;
      cout << "      - name with ID: '" << description2 << "'" << endl;
      cout << "      - name alone:   '" << description3 << "'" << endl;
      assert(description1 == description2);
      assert(description1 == description3);
    }

    cout << endl
         << "- " << _index << endl
         << string(80, '*') << endl
         << endl;

  }

  set<size_t> testIndexFile(const string &fname, DNAFileReader::Format fmt, DNAFileIndex::Status expected_status = DNAFileIndex::SUCCESS) {
    setFilename(fname, fmt);
    const string &new_fname = _reader.getFilename();
    DNAFileReader::Format fmt2 = _reader.getFormat();
    size_t expected_new_nb_nucleotides = 0;
    size_t expected_new_nb_sequences = 0;
    cout << "* Scanning the whole file (without indexing)..." << endl;
    set<size_t> bad_ids;
    list<string> indexed_sequences = _index.getSequenceNames(new_fname);
    size_t cpt = 0;
    bool check_consistency = _reader.check_consistency;
    bool warn = _reader.warn;
    assert(_index.consistencyChecking() == check_consistency);
    do {
      assert(_index.consistencyChecking() == check_consistency);
      assert(_reader.check_consistency == check_consistency);
      try {
        if (_reader.gotoNextSequence()) {
          size_t id = _reader.getCurrentSequenceID();
          ++cpt;
          assert(cpt == id);
          cout << "New sequence with ID " << id << " found." << endl;
          _reader.warn = false;
          _reader.computeCurrentSequenceLength();
          _reader.warn = warn;
          string name = DNAFileIndex::getNameFromDescription(_reader.getCurrentSequenceDescription());
          if (bad_ids.empty()) {
            expected_new_nb_nucleotides += _reader.getCurrentSequenceLength();
            if (find(indexed_sequences.begin(), indexed_sequences.end(), name) == indexed_sequences.end()) {
              ++expected_new_nb_sequences;
            }
          }
        }
      } catch (const FileReaderParseError &) {
        size_t id = _reader.getCurrentSequenceID();
        assert(cpt == id);
        bad_ids.insert(id);
        cout << "=> A parse error exception was caught for sequence " << id << ". <=" << endl;
        _reader.check_consistency = false;
        _reader.gotoSequenceEnd();
        _reader.check_consistency = check_consistency;
      }
    } while (_reader);
    _reader.close();

    cout << endl << "* Indexing the whole file..." << endl;
    DNAFileIndex::Status status;
    try {
      status = _index.indexFile(new_fname);
    } catch (const FileReaderParseError &) {
      size_t bad_id = _index.getNumberOfSequences(new_fname);
      assert(bad_ids.find(bad_id) != bad_ids.end());
      cout << "=> Expected FileReaderParseError caught for sequence " << bad_id << ". <=" << endl;
    }
    cout << "=> " << status << " (expecting " << expected_status << ")" << endl;
    fmt2 = _index.getFileFormat(new_fname);
    cout << "Now, detected format is '" << fmt2 << "'"
         << " (expecting '" << fmt << "')" << endl
         << endl;
    if (status == DNAFileIndex::SUCCESS) {
      _expected_nb_nucleotides += expected_new_nb_nucleotides;
      _expected_nb_sequences += expected_new_nb_sequences;
    } else {
      _expected_nb_nucleotides = _index.getNumberOfNucleotides();
      _expected_nb_sequences = _index.getNumberOfSequences();
    }
    assert(fmt == fmt2);
    return bad_ids;
  }

  void testIndexCurrentSeq(DNAFileIndex::Status expected_status = DNAFileIndex::SUCCESS) {
    cout << "* Indexing the current sequence" << endl;
    DNAFileIndex::Status status = _index.indexCurrentSequence(_reader);
    cout << "=> " << status << " (expecting " << expected_status << ")" << endl;
    assert(status == expected_status);

    _expected_nb_sequences += (status == DNAFileIndex::SUCCESS);

    // All sequence of the test file have the same pattern:
    regex pattern("[ ]*Sequence_([0-9]+) of length ([0-9]+) \\[([0-9]+) degenerated symbols\\] \\(from line ([0-9]+) to line ([0-9]+)\\)(.*)");
    smatch matches;
    // If NDEBUG macro is defined, the following code will not
    // compile, but this doesn't matter since the NDEBUG is undefined
    // at the beginning of this file.
    assert(regex_match(_reader.getCurrentSequenceDescription(), matches, pattern));
    assert(matches.size() == 7);
    size_t v = stoi(matches[1]);
    cout << "- Sequence ID (from the description string) is " << v << endl;
    cout << "- Sequence ID (from the reader) is " << _reader.getCurrentSequenceID() << endl;
    assert(v == _reader.getCurrentSequenceID());
    if (status == DNAFileIndex::SUCCESS) {
      size_t nb_nucl = stoi(matches[2]);
      cout << "- Sequence length is " << nb_nucl << endl;
      _expected_nb_nucleotides += nb_nucl;
    } else {
      cout << "Skipping the current sequence..." << endl;
      gotoNextSequence();
    }
    cout << endl;
    check();
  }

  void gotoNextSequence() {
    _reader.gotoNextSequence();
    _index.sync(_reader);
  }

  void testSet(const string &sequence_name, size_t pos, size_t expected_id,
                DNAFileIndex::Status expected_status = DNAFileIndex::SUCCESS,
                size_t nb_seq = 0, size_t nb_nucl = 0) {
    cout << "Moving reader to sequence '" << sequence_name << "' "
         << "at position " << pos << " (starting from 1)." << endl;
    assert(pos > 0);
    DNAFileIndex::Status status = _index.set(_reader, _reader.getFilename(), sequence_name, pos - 1);
    cout << "=> " << status << " (expecting " << expected_status << ")" << endl;
    assert(status == expected_status);
    _expected_nb_sequences += nb_seq;
    _expected_nb_nucleotides += nb_nucl;
    size_t id = _reader.getCurrentSequenceID();
    cout << "Sequence ID (from reader) is " << id
         << " (expecting " << expected_id << ")" << endl;
    assert(id == expected_id);

    check();
  }

  DNAFileReader::FileState getSequenceFileState() {
    return _reader.getState();
  }

  DNAFileReader::FileState getSequenceStartFileState() {
    DNAFileReader::FileState _expected_state = _reader.getState();
    _reader.gotoSequenceStart();
    DNAFileReader::FileState saved_state = _reader.getState();
    _reader.setState(_expected_state);
    assert(_reader.getState() == _expected_state);
    return saved_state;
  }

  const DNAFileReader &getReader() const {
    return _reader;
  }

  const DNAFileIndex &getIndex() const {
    return _index;
  }

};


int main() {

  DNAFileIndex_Tester tester;

  tester.consistencyChecking(true);
  tester.warn(true);

  // Test empty index.
  tester.check();

  // Index whole file to compute nb of nucl.
  set<size_t> bad_ids = tester.testIndexFile("test-reads.fastq", DNAFileReader::FASTQ_FORMAT);
  size_t nb_nucleotides = tester.getIndex().getNumberOfNucleotides();

  tester.check();

  // Clearing index.
  tester.clear();

  // Test empty index.
  tester.check();

  tester.setFilename("test-reads.fastq", DNAFileReader::FASTQ_FORMAT);

  tester.check();


  DNAFileReader::FileState saved_state;

  for (size_t i = 1; i <= 18; ++i) {
    try {
      assert(((tester.getReader().getCurrentSequenceID() == 0) && (i == 1))
             || ((i > 1) && (i == tester.getReader().getCurrentSequenceID())));
      tester.testIndexCurrentSeq();
      if (i == 1) {
        saved_state = tester.getSequenceStartFileState();
      }
      tester.gotoNextSequence();
    } catch (const FileReaderParseError &) {
      assert(bad_ids.find(i) != bad_ids.end());
      cout << "=> A parse error exception was caught for sequence " << i << ". <=" << endl;
      bad_ids.erase(i);
      // Test if nb of nucl. corresponds to the information obtained when
      // indexing the whole file.
      cout << "Number of nucleotides until this exception is " << tester.getIndex().getNumberOfNucleotides()
           << " (expecting " << nb_nucleotides << ")" << endl;
      assert(tester.getIndex().getNumberOfNucleotides() == nb_nucleotides);
      tester.consistencyChecking(false);
      tester.gotoNextSequence();
      tester.consistencyChecking(true);
    }
  }
  assert(bad_ids.empty());

  // Going to Sequence_1 at position 1 (starting from pos 1)
  tester.testSet("Sequence_1", 1, 1, DNAFileIndex::SUCCESS);

  DNAFileReader::FileState current_state = tester.getSequenceFileState();

  cout << "Current reader state is " << current_state << endl;
  cout << "Expecting to be in state " << saved_state << endl;
  assert(current_state == saved_state);

  // Going to Sequence_5 at position 20 (starting from pos 1)
  tester.testSet("Sequence_4", 20, 4, DNAFileIndex::SUCCESS);

  // Going to Sequence_1 at position 1 (starting from pos 1)
  tester.testSet("Sequence_1", 1, 1, DNAFileIndex::SUCCESS);

  current_state = tester.getSequenceFileState();

  cout << "Current reader state is " << current_state << endl;
  cout << "Expecting to be in state " << saved_state << endl;
  assert(current_state == saved_state);

  // Going to Sequence_4 at position 500 (starting from pos 1) which
  // doesn't exist
  tester.testSet("Sequence_4", 500, 4, DNAFileIndex::POSITION_NOT_FOUND);

  // Going to Sequence_1 at position 1 (starting from pos 1)
  tester.testSet("Sequence_1", 1, 1, DNAFileIndex::SUCCESS);

  current_state = tester.getSequenceFileState();

  cout << "Current reader state is " << current_state << endl;
  cout << "Expecting to be in state " << saved_state << endl;
  assert(current_state == saved_state);

  // Going to inexistant Sequence_1000 at position 50 (starting from
  // pos 1)
  tester.testSet("Sequence_1000", 50, 0, DNAFileIndex::SEQUENCE_NOT_FOUND);

  // Going to Sequence_1 at position 1 (starting from pos 1)
  tester.testSet("Sequence_1", 1, 1, DNAFileIndex::SUCCESS);

  current_state = tester.getSequenceFileState();

  cout << "Current reader state is " << current_state << endl;
  cout << "Expecting to be in state " << saved_state << endl;
  assert(current_state == saved_state);


  // Going to Sequence_10 at position 15 (starting from pos 1) and
  // saving file state
  tester.testSet("Sequence_10", 15, 10, DNAFileIndex::SUCCESS);
  saved_state = tester.getSequenceFileState();
  size_t nb_nucl_at_seq10 = 0;
  bool sequence_found = false;
  for (const auto &f: tester.getIndex().getFilenames()) {
    for (const auto &s: tester.getIndex().getSequenceNames(f)) {
      if (s == "Sequence_10") {
        sequence_found = true;
        break;
      }
      nb_nucl_at_seq10 += tester.getIndex().getNumberOfNucleotides(f, s);
    }
  }
  size_t nb_nucl_at_seq10_15 = nb_nucl_at_seq10 + (15 / tester.getIndex().bookmark_distance) * tester.getIndex().bookmark_distance;
  assert(sequence_found);
  cout << "Until 'Sequence_10' at position 15 there are " << nb_nucl_at_seq10_15 << " nucleotides." << endl;

  tester.clear();

  // test-sequences.fasta contains 6 sequences from Sequence_1 to
  // Sequence_6.
  tester.testIndexFile("test-sequences.fasta", DNAFileReader::FASTA_FORMAT);

  tester.check();

  tester.setFilename("test-reads.fastq", DNAFileReader::FASTQ_FORMAT);
  for (size_t i = 1; i <= 6; ++i) {
    // Sequence_1 to Sequence_6
    tester.testIndexCurrentSeq(DNAFileIndex::SEQUENCE_NAME_DUPLICATED);
  }

  tester.clear();
  tester.setFilename("test-reads.fastq", DNAFileReader::FASTQ_FORMAT);

  // Going to Sequence_10 at position 15 (starting from pos 1)
  tester.testSet("Sequence_10", 15, 10, DNAFileIndex::SUCCESS, 10, nb_nucl_at_seq10_15);

  current_state = tester.getSequenceFileState();

  // Since the fastq file (at least the 'Sequence_10') is not entirely
  // processed, the total number of nucleotides for this sequence will
  // differ. Thus forcing the correct total length of the sequence
  // prevents the assertion from failing for a bad reason.
  current_state.current_sequence_length = saved_state.current_sequence_length;
  cout << "Current reader state is " << current_state << endl;
  cout << "Expecting to be in state " << saved_state << endl;
  assert(current_state == saved_state);

  tester.testIndexFile("test-sequences.fasta", DNAFileReader::FASTA_FORMAT, DNAFileIndex::SEQUENCE_NAME_DUPLICATED);

  tester.setFilename("test-sequences.fasta", DNAFileReader::FASTA_FORMAT);

  for (size_t i = 1; i <= 6; ++i) {
    // Sequence_1 to Sequence_6
    tester.testIndexCurrentSeq(DNAFileIndex::SEQUENCE_NAME_DUPLICATED);
  }

  cout << "That's All, Folk!!!" << endl;
  return 0;

}
