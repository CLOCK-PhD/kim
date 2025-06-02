/******************************************************************************
*                                                                             *
*  Copyright © 2023-2025 -- IGH / LIRMM / CNRS / UM                           *
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
#include <dna_file_reader.h>

#include <cstdlib>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <regex>
#ifdef NDEBUG
#  undef NDEBUG
#endif
#include <cassert>

using namespace std;
using namespace kim;

void printSequenceInformations(const DNAFileReader &reader) {
  cout << "[" << __FUNCTION__ << "] "
       << "New sequence (with ID " << reader.getCurrentSequenceID() << "): "
       << "'" << reader.getCurrentSequenceDescription() << "'" << endl;
}

void printReaderDetails(const DNAFileReader &reader) {
  cout << "[" << __FUNCTION__ << "] "
       << reader.getFormat() << " file: '" << reader.getFilename() << "' ("
       << "line " << reader.getFileLineNumber() << ", "
       << "column " << reader.getFileColumnNumber() << ")."
       << endl;
}

void cancelAndExit(const DNAFileReader &) {
  cout << "[" << __FUNCTION__ << "] "
       << "This should not be invoked"
       << endl;
  exit(1);
}

// Check for opening & finding file methods
void test_DNAFileReader_init(string fname, DNAFileReader &reader, DNAFileReader::Format expected_fmt) {

  cout << "=== " << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__
       << "("
       << "fname = '" << fname << "', "
       << "reader@" << &reader << ", "
       << "check_consistency = " << reader.check_consistency << ", "
       << "expected_fmt = '" << expected_fmt << "'"
       << ") ==="
       << endl;

  reader.warn = true; // Turn on warnings

  cout << "Trying to open file '" << fname << "' which is not in the current working directory." << endl;
  reader.open(fname);
  cout << "Current opened filename should be empty: '" << reader.getFilename() << "'" << endl;
  assert(reader.getFilename().empty());
  cout << "Current file format should be Undefined: '" << reader.getFormat() << "'" << endl;
  assert(reader.getFormat() == DNAFileReader::UNDEFINED_FORMAT);

  cout << "Looking for file '" << fname << "' in the package directories." << endl;
  const char *dirs[] = {
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
  fname = FileReader::findFile(fname, dirs);
  cout << "File path should not be empty: '" << fname << "'" << endl;
  assert(!fname.empty());

  reader.open(fname);

  cout << "The current filename is '" << reader.getFilename() << "'" << endl;
  assert(reader.getFilename() == fname);
  cout << "Current file format should be " << expected_fmt << ": '" << reader.getFormat() << "'" << endl;
  assert(reader.getFormat() == expected_fmt);

  cout << "Closing file" << endl;
  reader.close();
  assert(!reader);

  cout << "Opening file again" << endl;
  reader.open(fname);
  assert(reader);

  cout << endl;
}

// test moving sequence by sequence in the file
vector<size_t> test_DNAFileReader_sequence_analysis(DNAFileReader &reader, size_t expected_nb_seq, list<size_t> bad_sequences) {

  cout << "=== " << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__
       << "("
       << "reader@" << &reader << ", "
       << "check_consistency = " << reader.check_consistency << ", "
       << "expected_nb_seq = " << expected_nb_seq << ", "
       << "bad_sequences = {";
  for (list<size_t>::const_iterator it = bad_sequences.begin();
       it != bad_sequences.end();
       ++it) {
    cout << (it == bad_sequences.begin() ? "" : ", ") << *it;
  }
  cout << "}"
       << ") ==="
       << endl;

  reader.warn = true; // Turn on warnings

  cout << "Resetting the reader before counting the number of sequences" << endl;
  reader.reset();
  assert(reader);

  vector<size_t> sequence_length;
  sequence_length.reserve(expected_nb_seq);
  size_t cpt = 0;

  if (reader.check_consistency) {

    // Check for proper but slower sequence counting This also check for
    // determining the sequence length and the save and restore cursor
    // positions.
    cout << "Counting the number of sequences checking for validity" << endl;
    while (reader) {

      try {

        ++cpt;
        if (!bad_sequences.empty() && (cpt == bad_sequences.front())) {
          assert(reader.getFormat() == DNAFileReader::FASTQ_FORMAT);
          cout << "The sequence " << bad_sequences.front()  << " has a longer quality and should throw an exception." << endl;
        }

        if (reader.gotoNextSequence()) {

          // All sequence of the test file have the same pattern:
          regex pattern("[ ]*Sequence_([0-9]+) of length ([0-9]+) \\[([0-9]+) degenerated symbols\\] \\(from line ([0-9]+) to line ([0-9]+)\\)(.*)");
          smatch matches;
          // If NDEBUG macro is defined, the following code will not
          // compile, but this doesn't matter since the NDEBUG is
          // undefined at the beginning of this file.
          assert(regex_match(reader.getCurrentSequenceDescription(), matches, pattern));
          assert(matches.size() == 7);
          size_t v = stoi(matches[1]);
          cout << "- Sequence ID (from the description string) is " << v << endl;
          cout << "- Sequence ID (from the reader) is " << reader.getCurrentSequenceID() << endl;
          assert(v == reader.getCurrentSequenceID());
          cout << "- Sequence ID (expected) is " << cpt << endl;
          assert(v == cpt);
          v = stoi(matches[2]);
          cout << "- Sequence length is " << v << endl;
          sequence_length.push_back(v);
          v = stoi(matches[4]);
          cout << "- Starting line is " << v << " (expecting " << reader.getFileLineNumber() << ")" << endl;
          assert(size_t(v) == reader.getFileLineNumber());

          DNAFileReader::FileState s1 = reader.getState();
          cout << "Before computing the sequence length, it should be not known... (expecting " << (ssize_t) reader.getCurrentSequenceLength() << ")" << endl;
          cout << "Sequence cursor is at " << reader.getFilename() << ":" << reader.getFileLineNumber() << ":" << reader.getFileColumnNumber() << endl;
          cout << "Start sequence cursor is at " << reader.getFilename() << ":" << (s1.sequence_start_file_state.line + 1) << ":" << (s1.sequence_start_file_state.column + 1) << "(" << s1.sequence_start_file_state.pos << ")" << endl;
          assert(reader.getCurrentSequenceLength() == size_t(-1));
          size_t length = reader.computeCurrentSequenceLength();
          cout << "After computation, the length of the sequence is " << length << endl;
          assert(length == sequence_length.back());
          cout << "Sequence cursor is now at " << reader.getFilename() << ":" << reader.getFileLineNumber() << ":" << reader.getFileColumnNumber() << endl;
          DNAFileReader::FileState s2 = reader.getState();
          assert(s1.pos == s2.pos);
          assert(s1.line == s2.line);
          assert(s1.column == s2.column);
          assert(s1.line + 1 == reader.getFileLineNumber());
          assert(s1.column + 1 == reader.getFileColumnNumber());
          assert(s1.start_symbol_expected == s2.start_symbol_expected);
          assert(s1.nb_nucleotides == s2.nb_nucleotides);
          assert(s1.nb_consecutive_regular_nucleotides == s2.nb_consecutive_regular_nucleotides);
          assert(s1.current_sequence_id == s2.current_sequence_id);
          assert(s1.current_sequence_description == s2.current_sequence_description);
          assert(s1.current_sequence_length != s2.current_sequence_length);
          assert(s1.kmer_aux == s2.kmer_aux);
          assert(s1.kmer_start_pos == s2.kmer_start_pos);
          assert(s1.kmer == s2.kmer);
          assert(s1.sequence_start_file_state.pos == s2.sequence_start_file_state.pos);
          assert(s1.sequence_start_file_state.line == s2.sequence_start_file_state.line);
          assert(s1.sequence_start_file_state.column == s2.sequence_start_file_state.column);
          assert(s1.sequence_start_file_state.pos == s2.pos);
          assert(s1.sequence_start_file_state.line == s2.line);
          assert(s1.sequence_start_file_state.column == s2.column);
        } else {
          cout << "No new sequence found." << endl;
          --cpt;
        }

      } catch (const FileReaderParseError &e) {

        cout << "A FileReaderParseError with message '" << e.what() << "' has been thrown" << endl;
        if (!bad_sequences.empty() && (cpt == bad_sequences.front() + 1)) {
          assert(reader.getFormat() == DNAFileReader::FASTQ_FORMAT);
          bad_sequences.pop_front();
          cout << "This is expected since this sequence is voluntary badly formatted." << endl;
          cout << "Skip the end of this quality sequence until a new sequence start is being detected..." << endl;
          reader.check_consistency = false;
          reader.gotoSequenceEnd();
          reader.check_consistency = true;
          --cpt;
        } else {
          exit(1);
        }

      }

    }

    assert(bad_sequences.empty());

  } else {

    // Check for quick and dirty sequence counting
    cout << "Counting the number of sequences using no checking (thus will be wrong)" << endl;
    while (reader && (cpt < 1000)) {
      if (reader.gotoNextSequence()) {
        ++cpt;
        cout << "- Sequence ID (from the reader) is " << reader.getCurrentSequenceID() << endl;
        cout << "- Sequence ID (expected) is " << cpt << endl;
        assert(cpt == reader.getCurrentSequenceID());
      } else {
        cout << "No new sequence found." << endl;
      }
    }

  }

  cout << "There is " << cpt << " lines (expecting " << expected_nb_seq << ") seen as sequences in '" << reader.getFilename() << "'" << endl;
  assert(cpt == expected_nb_seq);

  cout << endl;

  return sequence_length;

}


// Check for k-mers extraction and analysis methods.
void test_DNAFileReader_kmer_extraction(DNAFileReader &reader, bool skip_degenerated_kmers, const vector<size_t> &sequence_length, list<size_t> bad_sequences, vector<list<size_t>> degenerated_kmers) {
  cout << "=== " << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__
       << "("
       << "reader@" << &reader << ", "
       << "check_consistency = " << reader.check_consistency << ", "
       << "skip_degenerated_kmers = " << skip_degenerated_kmers << ", "
       << "sequence_length = [";
  bool check_consistency = reader.check_consistency;
  for (vector<size_t>::const_iterator it = sequence_length.begin();
       it != sequence_length.end();
       ++it) {
    cout << (it == sequence_length.begin() ? "" : ", ") << *it;
  }
  cout << "], "
       << "bad_sequences = {";
  for (list<size_t>::const_iterator it = bad_sequences.begin();
       it != bad_sequences.end();
       ++it) {
    cout << (it == bad_sequences.begin() ? "" : ", ") << *it;
  }
  cout << "}, "
       << "degenerated_kmers = [";
  for (size_t i = 1; i < degenerated_kmers.size(); ++i) {
    cout << (i == 1 ? "" : ", ") << "{";
    for (list<size_t>::const_iterator it = degenerated_kmers[i].begin();
         it != degenerated_kmers[i].end();
         ++it) {
      cout << (it == degenerated_kmers[i].begin() ? "" : ", ") << *it;
    }
    cout << "}";
  }
  cout << "]"
       << ") ==="
       << endl;

  reader.warn = false; // Turn off warnings

  cout << "Resetting the reader before analyzing k-mers." << endl;
  reader.reset();
  assert(reader);

  size_t last_id = 0;
  size_t k = reader.k();
  size_t nb_kmers = -1;
  for (string kmer = reader.getNextKmer(skip_degenerated_kmers); !kmer.empty(); kmer = reader.getNextKmer(skip_degenerated_kmers)) {

    size_t p = reader.getCurrentKmerRelativePosition();
    size_t cpt = reader.getCurrentSequenceID();
    if (last_id != cpt) {
      if (last_id) {
        cout << "End of the sequence " << last_id << " of length " << sequence_length[last_id - 1]
             << " (" << nb_kmers << " processed or skipped k-mers)" << endl;
        assert(nb_kmers + k - 1 == sequence_length[last_id - 1]);
      } else {
        assert(nb_kmers == size_t(-1));
      }
      if (skip_degenerated_kmers) {
        while (++last_id < cpt) {
          degenerated_kmers[last_id].clear();
        }
      } else {
        last_id = cpt;
      }
      if (degenerated_kmers[cpt].empty()) {
        cout << "All k-mers from this sequence contains no degenerated symbols" << endl;
      } else {
        cout << "List of k-mers that contains degenerated symbols in sequence " << cpt << ":" << endl;
        for (list<size_t>::const_iterator it = degenerated_kmers[cpt].begin();
             it != degenerated_kmers[cpt].end();
             ++it) {
          cout << (it == degenerated_kmers[cpt].begin() ? "" : ", ") << *it;
        }
      }
      nb_kmers = 0;
    }
    cout << "Current " << k << "-mer '" << kmer << "'  is at position " << p << endl;
    cout << "(reader.getCurrentKmer() = '" << reader.getCurrentKmer() << "' and sequence length is " << sequence_length[cpt - 1] << ")" << endl;
    assert(kmer == reader.getCurrentKmer());
    assert(p + k <= sequence_length[cpt - 1]);
    cout << "Checking if current k-mer contains degenerated symbols: "
         << reader.currentKmerContainsDegenerateSymbols() << endl;
    if (!degenerated_kmers[cpt].empty() && !skip_degenerated_kmers) {
      if (p == degenerated_kmers[cpt].front()) {
        assert(reader.currentKmerContainsDegenerateSymbols());
        degenerated_kmers[cpt].pop_front();
      } else {
        assert(!reader.currentKmerContainsDegenerateSymbols());
      }
    } else {
      while (!degenerated_kmers[cpt].empty() && (degenerated_kmers[cpt].front() < p)) {
        degenerated_kmers[cpt].pop_front();
        ++nb_kmers;
      }
      assert(degenerated_kmers[cpt].empty() || (degenerated_kmers[cpt].front() > p));
      assert(!reader.currentKmerContainsDegenerateSymbols());
    }
    cout << "Expected k-mer position is " << nb_kmers << endl;
    assert(p == nb_kmers);
    ++nb_kmers;

    if (!bad_sequences.empty() && (cpt == bad_sequences.front()) && (nb_kmers + k == sequence_length[cpt - 1])) {
      cout << "In order to not throw and exception, the method gotoSequenceEnd() is invoked." << endl;
      reader.check_consistency = false;
      reader.gotoSequenceEnd();
      reader.check_consistency = check_consistency;
      bad_sequences.pop_front();
      ++nb_kmers;
    }
  }
  cout << "End of file after having processed " << reader.getCurrentSequenceID() << " sequences (expecting " << sequence_length.size() << ")" << endl;
  assert(reader.getCurrentSequenceID() == sequence_length.size());

  assert(bad_sequences.empty());
  for (size_t i = 0; i < degenerated_kmers.size(); ++i) {
    cout << "Ensure that all degenerated kmers have been encountered for sequence " << i << endl;
    cout << "degenerated_kmers[" << i << "] = {";
    for (list<size_t>::const_iterator it = degenerated_kmers[i].begin();
         it != degenerated_kmers[i].end();
         ++it) {
      cout << (it == degenerated_kmers[i].begin() ? "" : ", ") << *it;
    }
    cout << "}" << endl;
    assert(degenerated_kmers[i].empty());
  }

  cout << endl;

}

void test_DNAFileReader_indexation(DNAFileReader &reader, list<size_t> bad_sequences) {
  cout << "=== " << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__
       << "("
       << "reader@" << &reader << ", "
       << "check_consistency = " << reader.check_consistency << ", "
       << "bad_sequences = {";
  for (list<size_t>::const_iterator it = bad_sequences.begin();
       it != bad_sequences.end();
       ++it) {
    cout << (it == bad_sequences.begin() ? "" : ", ") << *it;
  }
  cout << "}"
       << ") ==="
       << endl;

  reader.warn = false; // Turn off warnings

  cout << "Resetting the reader before building a quick sequence index." << endl;
  reader.reset();
  assert(reader);

  vector<list<DNAFileReader::FileState>> index;
  vector<list<string>> kmers;

  cout << "  Building the sequence index" << endl;
  size_t last_id = reader.getCurrentSequenceID();
  size_t p = 0;
  while (reader) {
    try {
      if (reader.gotoNextSequence()) {
        if (last_id != reader.getCurrentSequenceID()) {
          last_id = reader.getCurrentSequenceID();
          cout << "    Indexation of sequence " << last_id << endl;
          index.emplace_back();
          kmers.emplace_back();
          cout << "    Starting a new list at index[" << index.size() - 1 << "]" << endl;
          assert(index.size() == last_id);
          p = 0;
        }

        list<DNAFileReader::FileState> &cur_states = index.back();
        list<string> &cur_kmers = kmers.back();
        while (reader.getCurrentSequenceLength() == size_t(-1)) {
          // The reader has not encountered the end of the current sequence.
          string kmer = reader.getKmerAt(p);
          if (reader.getCurrentSequenceLength() == (size_t)-1) {
            DNAFileReader::FileState state = reader.getState();
            cout << "      kmer at " << reader.getCurrentKmerRelativePosition()
                 << " (wanted position was " << p << ")"
                 << " corresponds to kmer '" << kmer
                 << "' from '" << state.filename
                 << ":" << state.line
                 << ":" << state.column
                 << "@" << state.pos << endl;
            assert(p == reader.getCurrentKmerRelativePosition());
            cur_kmers.push_back(kmer);
            cur_states.push_back(state);
            p += 25; // index every 25 positions
          }
        }
      }
    } catch (const FileReaderParseError &e) {

      cout << "A FileReaderParseError with message '" << e.what() << "' has been thrown" << endl;
      if (!bad_sequences.empty() && (last_id == bad_sequences.front())) {
        assert(reader.getFormat() == DNAFileReader::FASTQ_FORMAT);
        bad_sequences.pop_front();
        cout << "This is expected since this sequence is voluntary badly formatted." << endl;
        cout << "Skip the end of this quality sequence until a new sequence start is being detected..." << endl;
        bool check_consistency = reader.check_consistency;
        reader.check_consistency = false;
        reader.gotoSequenceEnd();
        reader.check_consistency = check_consistency;;
      } else {
        exit(1);
      }
    }
  }

  cout << "Validation of the sequence index" << endl;
  assert(index.size() == kmers.size());
  for (size_t i = 0; i < index.size(); ++i) {
    cout << "  Index of sequence " << (i + 1) << endl;
    list<DNAFileReader::FileState> &cur_states = index[i];
    list<string> &cur_kmers = kmers[i];
    assert(cur_states.size() == cur_kmers.size());
    while (!cur_kmers.empty()) {
      string expected_kmer = cur_kmers.front();
      cur_kmers.pop_front();
      DNAFileReader::FileState expected_state = cur_states.front();
      cur_states.pop_front();
      reader.setState(expected_state);
      DNAFileReader::FileState state = reader.getState();
      string kmer = reader.getCurrentKmer();
      size_t pos = reader.getCurrentKmerRelativePosition();
      cout << "      kmer at " << pos
           << " corresponds to kmer '" << kmer << "'"
           << " (expecting '" << expected_kmer << "')"
           << " from '" << state.filename
           << ":" << state.line
           << ":" << state.column
           << "@" << state.pos
           << " (expecting '" << expected_state.filename
           << ":" << expected_state.line
           << ":" << expected_state.column
           << "@" << expected_state.pos
           << ")" << endl;
      assert(kmer == expected_kmer);
      assert(state.filename == expected_state.filename);
      assert(state.line == expected_state.line);
      assert(state.column == expected_state.column);
      assert(state.pos == expected_state.pos);
    }
  }

  cout << endl;

}

struct kmer_record {
  string kmer;
  size_t pos;
  kmer_record(const string &s, size_t p):
    kmer(s), pos(p)
  {}
};

void test_DNAFileReader_some_kmers(DNAFileReader &reader, size_t num_seq, list<kmer_record> kmers) {

  cout << "=== " << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__
       << "("
       << "reader@" << &reader << ", "
       << "check_consistency = " << reader.check_consistency << ", "
       << "num_seq = " << num_seq << ", "
       << "kmers = {";
  for (list<kmer_record>::const_iterator it = kmers.begin();
       it != kmers.end();
       ++it) {
    cout << (it == kmers.begin() ? "" : ", ")
         << "('" << it->kmer << "', " << it->pos << ")";
  }
  cout << "]"
       << ") ==="
       << endl;

  reader.warn = false; // Turn off warnings

  cout << "Resetting the reader before in-depth analyzis of some k-mers." << endl;
  reader.reset();
  assert(reader);
  size_t k = reader.k();

  for (size_t i = 0; i < num_seq; ++i) {
      reader.gotoNextSequence();
  }

  for (list<kmer_record>::iterator it = kmers.begin();
       it != kmers.end();
       ++it) {
    cout << "[BEFORE]:" << endl
         << "k-mer: '" << reader.getCurrentKmer() << "'" << endl
         << "k-mer pos: " << reader.getCurrentKmerRelativePosition() << endl
         << "Seq ID: " << reader.getCurrentSequenceID() << endl
         << "Seq Description: '" << reader.getCurrentSequenceDescription() << "'" << endl
         << "Seq nb nucl: " << reader.getCurrentSequenceProcessedNucleotides() << endl
         << "Seq length: " << reader.getCurrentSequenceLength() << endl;

    reader.getKmerAt(it->pos);

    cout << "[AFTER]:" << endl
         << "k-mer: '" << reader.getCurrentKmer() << "'" << endl
         << "k-mer pos: " << reader.getCurrentKmerRelativePosition() << endl
         << "Seq ID: " << reader.getCurrentSequenceID() << endl
         << "Seq Description: '" << reader.getCurrentSequenceDescription() << "'" << endl
         << "Seq nb nucl: " << reader.getCurrentSequenceProcessedNucleotides() << endl
         << "Seq length: " << reader.getCurrentSequenceLength() << endl;


    cout << "The " << k << "-mer at position " << it->pos << " in sequence " << num_seq << " is '" << reader.getCurrentKmer() << "' (expecting '" << it->kmer << "')." << endl;
    assert(reader.getCurrentKmer() == it->kmer);
    cout << "It is located at position " << reader.getCurrentKmerRelativePosition() << " (expecting " << (it->kmer.empty() ? -1 : it->pos) << ")" << endl;
    assert(reader.getCurrentKmerRelativePosition() == (it->kmer.empty() ? -1 : it->pos));
  }

  cout << endl;

}


void _fill_list(list<size_t> &l, size_t first, size_t last) {
  while (first != last) {
    l.push_back(first);
    ++first;
  }
}

int main() {

  Settings settings(5, 2);
  assert(settings.valid());
  settings.freeze();
  assert(settings.frozen());
  cout << "Testing DNAFileReader"
       << " with k = " << settings.getKmerLength()
       << " and p = " << settings.getKmerPrefixLength()
       << endl;

  DNAFileReader reader(settings.k());
  cout << "Add two callback functions called each time a new sequence is encountered" << endl;
  reader.addOnNewSequenceCallback(printSequenceInformations);
  reader.addOnNewSequenceCallback(printReaderDetails);
  reader.addOnNewSequenceCallback(cancelAndExit);
  reader.removeOnNewSequenceCallback(cancelAndExit);
  vector<size_t> sequence_length;
  list<size_t> bad_sequences;
  vector<list<size_t>> degenerated_kmers;

  test_DNAFileReader_init("test-reads.fastq", reader, DNAFileReader::FASTQ_FORMAT);
  bad_sequences.push_back(15);
  reader.check_consistency = false;
  test_DNAFileReader_sequence_analysis(reader, 21, bad_sequences);
  reader.check_consistency = true;
  sequence_length = test_DNAFileReader_sequence_analysis(reader, 18, bad_sequences);
  degenerated_kmers.resize(18 + 1); // degenerated k-mers for sequence i are stored at offset i

  _fill_list(degenerated_kmers[4], 5, 10);
  _fill_list(degenerated_kmers[4], 36, 50);
  _fill_list(degenerated_kmers[5], 15, 25);
  _fill_list(degenerated_kmers[5], 35, 43);
  _fill_list(degenerated_kmers[5], 55, 60);
  _fill_list(degenerated_kmers[17], 0, 26);
  test_DNAFileReader_kmer_extraction(reader, true /* skip degenerated k-mers */, sequence_length, bad_sequences, degenerated_kmers);
  test_DNAFileReader_kmer_extraction(reader, false /* don't skip degenerated k-mers */, sequence_length, bad_sequences, degenerated_kmers);
  test_DNAFileReader_indexation(reader, bad_sequences);

  test_DNAFileReader_init("test-reads.fastq.gz", reader, DNAFileReader::FASTQ_FORMAT);
  test_DNAFileReader_kmer_extraction(reader, true /* skip degenerated k-mers */, sequence_length, bad_sequences, degenerated_kmers);
  test_DNAFileReader_kmer_extraction(reader, false /* don't skip degenerated k-mers */, sequence_length, bad_sequences, degenerated_kmers);
  test_DNAFileReader_indexation(reader, bad_sequences);

  test_DNAFileReader_init("test-reads.fastq.bz2", reader, DNAFileReader::FASTQ_FORMAT);
  test_DNAFileReader_kmer_extraction(reader, true /* skip degenerated k-mers */, sequence_length, bad_sequences, degenerated_kmers);
  test_DNAFileReader_kmer_extraction(reader, false /* don't skip degenerated k-mers */, sequence_length, bad_sequences, degenerated_kmers);
  test_DNAFileReader_indexation(reader, bad_sequences);

  sequence_length.clear();
  bad_sequences.clear();
  degenerated_kmers.clear();

  test_DNAFileReader_init("test-sequences.fasta", reader, DNAFileReader::FASTA_FORMAT);
  reader.check_consistency = false;
  test_DNAFileReader_sequence_analysis(reader, 6, bad_sequences);
  reader.check_consistency = true;
  sequence_length = test_DNAFileReader_sequence_analysis(reader, 6, bad_sequences);
  degenerated_kmers.resize(6 + 1); // same as above
  _fill_list(degenerated_kmers[3], 0, 19);
  _fill_list(degenerated_kmers[3], 20, 25);
  _fill_list(degenerated_kmers[4], 16, 38);
  _fill_list(degenerated_kmers[4], 44, 50);
  test_DNAFileReader_kmer_extraction(reader, true /* skip degenerated k-mers */, sequence_length, bad_sequences, degenerated_kmers);
  test_DNAFileReader_kmer_extraction(reader, false /* don't skip degenerated k-mers */, sequence_length, bad_sequences, degenerated_kmers);
  test_DNAFileReader_indexation(reader, bad_sequences);

  list<kmer_record> some_kmers;
  // The current sequence is the one with TuX
  some_kmers.emplace_back("ACTAC", 10);
  some_kmers.emplace_back("ACTAC", 10);
  some_kmers.emplace_back("CTACT", 11);
  some_kmers.emplace_back("TCGAT", 38);
  some_kmers.emplace_back("ATCNN", 45);
  some_kmers.emplace_back("ACTAC", 10);
  some_kmers.emplace_back("", 1000);
  some_kmers.emplace_back("ACTAC", 10);
  test_DNAFileReader_some_kmers(reader, 4, some_kmers);

  cout << "That's All, Folk!!!" << endl;
  return 0;

}
