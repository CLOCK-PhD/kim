/******************************************************************************
*                                                                             *
*  Copyright © 2023-2024 -- IGH / LIRMM / CNRS / UM                           *
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

#include "dna_file_reader.h"

#include "config.h"
#include "kim_exception.h"

#include <cassert>
#include <iostream>
#include <limits>

using namespace std;

BEGIN_KIM_NAMESPACE

#define WARNING_MSG(msg)                        \
  if (warn) {                                   \
    cerr << "Warning:"                          \
         << getFilename() << ":"                \
         << getFileLineNumber() << ":"          \
         << getFileColumnNumber() << ": "       \
         << msg << endl;                        \
  }                                             \
  (void) 0

#define ERROR_MSG(msg)                          \
  do {                                          \
    FileReaderParseError error(*this);          \
    error << msg;                               \
    throw error;                                \
  } while (0)



enum IUPAC {
            IUPAC_UNDEFINED = 0, // Undefined
            IUPAC_A = 1, // Adenine
            IUPAC_C = 2, // Cytosine
            IUPAC_G = 4, // Guanine
            IUPAC_T = 8, // Thymine
            IUPAC_U = 8, // Uracil
            IUPAC_GAP = 16, // A gap symbol (no nucleotide)
            IUPAC_R = IUPAC_A | IUPAC_G, // Purines
            IUPAC_Y = IUPAC_C | IUPAC_T, // Pyrimidines
            IUPAC_K = IUPAC_G | IUPAC_T, // Ketones
            IUPAC_M = IUPAC_A | IUPAC_C, // Amino groups
            IUPAC_S = IUPAC_C | IUPAC_G, // Strong interactions
            IUPAC_W = IUPAC_A | IUPAC_T, // Weak interactions
            IUPAC_B = IUPAC_C | IUPAC_G | IUPAC_T, // All but Adenine
            IUPAC_D = IUPAC_A | IUPAC_G | IUPAC_T, // All but Cytosine
            IUPAC_H = IUPAC_A | IUPAC_C | IUPAC_T, // All but Guanine
            IUPAC_V = IUPAC_A | IUPAC_C | IUPAC_G, // All but Thymine nor Uracil
            IUPAC_N = IUPAC_A | IUPAC_C | IUPAC_G | IUPAC_T // Any
};

static IUPAC toIUPAC(char c) {
  if ((c >= 'a') && (c <= 'z')) {
    c += 'A' - 'a'; // Upcase
  }
  switch (c) {
  case 'A': return IUPAC_A;
  case 'C': return IUPAC_C;
  case 'G': return IUPAC_G;
  case 'T': return IUPAC_T;
  case 'U': return IUPAC_U;
  case 'R': return IUPAC_R;
  case 'Y': return IUPAC_Y;
  case 'K': return IUPAC_K;
  case 'M': return IUPAC_M;
  case 'S': return IUPAC_S;
  case 'W': return IUPAC_W;
  case 'B': return IUPAC_B;
  case 'D': return IUPAC_D;
  case 'H': return IUPAC_H;
  case 'V': return IUPAC_V;
  case 'N': return IUPAC_N;
  case '.':
  case '-': return IUPAC_GAP; // Gap is valid in IUPAC even if it's
                              // not a nucleotide
  default:
    return IUPAC_UNDEFINED;
  }
}


DNAFileReader::FileState::FileState():
  FileReader::FileState(),
  format(DNAFileReader::UNDEFINED_FORMAT),
  k(0)
{
  reset();
}

void DNAFileReader::FileState::reset() {
  start_symbol_expected = true;
  nb_nucleotides = nb_consecutive_regular_nucleotides = 0;
  current_sequence_id = 0;
  current_sequence_description = "";
  current_sequence_length = size_t(-1);
  kmer = kmer_aux = string(k, '?');
  kmer_start_pos = 0;
  sequence_start_file_state = FileReader::FileState();
}


// Template instanciations (must be instantiated before being used...)

template <DNAFileReader::Format>
bool _isStartSymbol(char c);

template <>
bool _isStartSymbol<DNAFileReader::FASTA_FORMAT>(char c) {
  return ((c == '>') || (c == ';'));
}

template <>
bool _isStartSymbol<DNAFileReader::FASTQ_FORMAT>(char c) {
  return (c == '@');
}

#ifndef __UNUSED__
#  define __UNUSED__(x)
#endif

template <DNAFileReader::Format>
bool _isStartComment(char __UNUSED__(c)) {
  return false;
}

template <>
bool _isStartComment<DNAFileReader::FASTA_FORMAT>(char c) {
  return (c == ';');
}

template<DNAFileReader::Format fmt>
string start_symbols;

template<>
string start_symbols<DNAFileReader::FASTA_FORMAT> = ">'";

template<>
string start_symbols<DNAFileReader::FASTQ_FORMAT> = "@";

template <DNAFileReader::Format>
bool _isEndOfNucleotideSequence(size_t col, char c);

template <>
bool _isEndOfNucleotideSequence<DNAFileReader::FASTA_FORMAT>(size_t col, char c) {
  return ((c == -1) || ((col == 1) && _isStartSymbol<DNAFileReader::FASTA_FORMAT>(c)));
}

template <>
bool _isEndOfNucleotideSequence<DNAFileReader::FASTQ_FORMAT>(size_t col, char c) {
  return ((col == 1) && (c == '+'));
}

template<DNAFileReader::Format fmt>
bool DNAFileReader::_parseSequenceDescription() {
  int c = _nextVisibleCharacter();
  if (c) {
    FileState &current_state = _getState();
    if (current_state.column == 1) {
      if (c != 0) {
        if (_isStartSymbol<fmt>(c)) {
          _processSequenceDescription();
          current_state.sequence_start_file_state = current_state;
          for (list<onNewSequenceFct>::const_iterator it = _on_new_sequence_cb.begin();
               it != _on_new_sequence_cb.end();
               ++it) {
            (*it)(*this);
          }
        } else {
          ERROR_MSG("Badly formatted file. Character '" << (char) c << "'" << " found while expecting "
                    << ((start_symbols<fmt>.size() == 1) ? "'" : "one of '")
                    << start_symbols<fmt> << "'.");
        }
      } else {
        ERROR_MSG("Badly formatted file. End of file found while expecting "
                  << ((start_symbols<fmt>.size() == 1) ? "'" : "one of '")
                  << start_symbols<fmt> << "'.");
      }
    } else {
      ERROR_MSG("Badly formatted file. Unexpected character '" << (char) c << "'.");
    }
  }
  return c; // (c != '\0') => Correct Header
}

template <>
void DNAFileReader::_parseEndOfNucleotideSequence<DNAFileReader::FASTA_FORMAT>() {
  FileState &current_state = _getState();
  if (_ifs) {
    _ifs.unget();
    assert(current_state.column > 0);
    --current_state.column;
  } else {
    _ifs.clear();
  }
  current_state.pos = _ifs.tellg();
  assert((current_state.current_sequence_length == size_t(-1)) ||
         (current_state.current_sequence_length == current_state.nb_nucleotides));
  current_state.current_sequence_length = current_state.nb_nucleotides;
}

template <>
void DNAFileReader::_parseEndOfNucleotideSequence<DNAFileReader::FASTQ_FORMAT>() {
  // Handle the separator line between nucleotides and quality
  // symbols.
  string description;
  getline(_ifs, description);
  FileState &current_state = _getState();
  ++current_state.line;
  current_state.column = 0;
  current_state.pos = _ifs.tellg();
  if (!description.empty() && (description != current_state.current_sequence_description)) {
    WARNING_MSG("Badly formatted file. The description following the '+' sign (" << description << ") is not the expected one (" << current_state.current_sequence_description << ").");
  }

  // Now gobble the quality
  int c = -1;
  size_t nb = current_state.nb_nucleotides;
  while (c && nb--) {
    c = _nextVisibleCharacter();
  }
  if (!c) {
    ERROR_MSG("Badly formatted file. End of file found while " << nb << " quality symbols were still expected.");
  }
  c = _nextVisibleCharacter();
  _parseEndOfNucleotideSequence<FASTA_FORMAT>();
}

template <DNAFileReader::Format fmt>
void DNAFileReader::_parse() {
  bool ok;
  FileState &current_state = _getState();

  if (current_state.start_symbol_expected) {
    ok = _parseSequenceDescription<fmt>();
  } else {
    do {
      int c = _nextVisibleCharacter();
      if (c) {
        IUPAC nucl = toIUPAC(c);
        size_t p = current_state.kmer_start_pos;
        if (current_state.nb_nucleotides < current_state.k) {
          p += current_state.nb_nucleotides;
        }
        // cerr << "The next valid symbol (" << (char) c << "?) will be added at position " << p << " since kmer starts at position " << current_state.kmer_start_pos << endl;
        switch (nucl) {
        case IUPAC_GAP:
          // ignore gaps
          ok = false;
          break;
        case IUPAC_UNDEFINED:
          if (_isEndOfNucleotideSequence<fmt>(current_state.column, c)) {
            _parseEndOfNucleotideSequence<fmt>();
            current_state.start_symbol_expected = true;
            ok = true;
          } else {
            if (_isStartComment<fmt>(c)) {
              _ifs.ignore(numeric_limits<streamsize>::max(), '\n');
              ++current_state.line;
              current_state.column = 0;
              current_state.pos = _ifs.tellg();
            } else {
              WARNING_MSG("Badly formatted file. Unexpected character '" << (char) c << "'.");
            }
          }
          break;
        default:
          // This is a valid symbol
          c = toupper(c);
          current_state.kmer_aux[p] = c;
          if (current_state.nb_nucleotides++ >= current_state.k) {
            if (++current_state.kmer_start_pos >= current_state.k) {
              current_state.kmer_start_pos = 0;
            }
          }
          // cerr << "Now sequence length is " << current_state.nb_nucleotides << endl;
          ok = (current_state.nb_nucleotides >= current_state.k);
          switch (nucl) {
          case IUPAC_A:
          case IUPAC_C:
          case IUPAC_G:
          case IUPAC_T:
            ++current_state.nb_consecutive_regular_nucleotides;
            break;
          default:
            // This is a degenerate symbol
            current_state.nb_consecutive_regular_nucleotides = 0;
          }
        }
      } else {
        assert(!*this);
        ok = false;
        assert((current_state.current_sequence_length == size_t(-1)) ||
               (current_state.current_sequence_length == current_state.nb_nucleotides));
        current_state.current_sequence_length = current_state.nb_nucleotides;
      }
    } while (*this && !ok);
  }
  assert(ok xor !*this);
}



DNAFileReader::DNAFileReader(size_t k, const string &filename, bool warn):
  FileReader("", warn)
{
  _getState().k = k;
  if (!filename.empty()){
    open(filename);
  }
}

void DNAFileReader::_onOpen() {
  FileState &current_state = _getState();
  current_state.format = UNDEFINED_FORMAT;
  _parse_mth = NULL;
  if (!_ifs) return;
  char c = _nextVisibleCharacter();
  if (_isStartSymbol<FASTA_FORMAT>(c)) {
    current_state.format = FASTA_FORMAT;
    _parse_mth = &DNAFileReader::_parse<FASTA_FORMAT>;
  } else if (_isStartSymbol<FASTQ_FORMAT>(c)) {
    current_state.format = FASTQ_FORMAT;
    _parse_mth = &DNAFileReader::_parse<FASTQ_FORMAT>;
  }
  if (current_state.format != UNDEFINED_FORMAT) {
    _ifs.unget();
    assert(current_state.column > 0);
    --current_state.column;
    current_state.pos = _ifs.tellg();
  } else {
    ERROR_MSG("Unable to detect the file format. Character '" << (char) c << "'"
              << " found while expecting either a '@' for fastq file or one of '>' or ';' for fasta file.");
  }
}

void DNAFileReader::_onClose() {
  _getState().format = UNDEFINED_FORMAT;
}

void DNAFileReader::_onReset() {
  FileState &current_state = _getState();
  current_state.reset();
  current_state.sequence_start_file_state = FileReader::getState();
}

void DNAFileReader::_processSequenceDescription() {
  FileState &current_state = _getState();
  current_state.start_symbol_expected = false;
  getline(_ifs, current_state.current_sequence_description);
  ++current_state.current_sequence_id;
  ++current_state.line;
  current_state.column = 0;
  current_state.pos = _ifs.tellg();
  current_state.nb_nucleotides = current_state.nb_consecutive_regular_nucleotides = 0;
  current_state.current_sequence_length = size_t(-1);
  current_state.kmer_start_pos = 0;
}

const string &DNAFileReader::getCurrentKmer() const {
  static const string empty_string;
  const FileState &current_state = getState();
  if (*this && (current_state.nb_nucleotides >= current_state.k)) {
    for (size_t p = current_state.kmer_start_pos; p < current_state.k; ++p) {
      current_state.kmer[p - current_state.kmer_start_pos] = current_state.kmer_aux[p];
    }
    for (size_t p = 0; p < current_state.kmer_start_pos; ++p) {
      current_state.kmer[p + current_state.k - current_state.kmer_start_pos] = current_state.kmer_aux[p];
    }
  } else {
    current_state.kmer = empty_string;
  }
  return current_state.kmer;
}

const string &DNAFileReader::getNextKmer(bool skip_degenerate) {
  size_t nb;
  const FileState &current_state = getState();
  do {
    (this->*_parse_mth)();
    if (current_state.start_symbol_expected) {
      (this->*_parse_mth)();
    }
    nb = (skip_degenerate ? current_state.nb_consecutive_regular_nucleotides : current_state.nb_nucleotides);
  } while (*this && (nb < current_state.k));
  return getCurrentKmer();
}

const string &DNAFileReader::getKmerAt(size_t p, bool check_consistency) {
  if (!*this) return getCurrentKmer();
  const FileState &current_state = getState();
  if (p < getCurrentKmerRelativePosition()) {
    gotoSequenceStart();
  }
  if (getCurrentKmerRelativePosition() == (size_t) -1) {
    assert(current_state.nb_nucleotides == 0);
    assert(current_state.nb_consecutive_regular_nucleotides == 0);
    // The sequence description was already processed, but no k-mer is
    // currently loaded.
    (this->*_parse_mth)();
  }
  assert((p >= getCurrentKmerRelativePosition())
         || (current_state.current_sequence_length < current_state.k));
  if (check_consistency || (p < current_state.nb_nucleotides)) {
    while (*this && !current_state.start_symbol_expected && (getCurrentKmerRelativePosition() < p)) {
      (this->*_parse_mth)();
    }
  } else {
    assert(p >= current_state.nb_nucleotides);
    while (*this && ((current_state.nb_nucleotides + 1) < p)) {
      _nextVisibleCharacter();
    }
    if (*this) {
      assert((current_state.nb_consecutive_regular_nucleotides + 1) == p);
      (this->*_parse_mth)();
    }
  }
  return getCurrentKmer();
}

template <DNAFileReader::Format fmt>
void DNAFileReader::_gotoSequenceEnd(bool check_consistency) {
  FileState &current_state = _getState();
  if (check_consistency) {
    while (*this && !current_state.start_symbol_expected) {
      _parse<fmt>();
    }
  } else {
    bool found = false;
    while (*this && !found) {
      int c = _nextVisibleCharacter();
      if ((current_state.column == 1) and _isStartSymbol<fmt>(c)) {
        found = true;
      }
    }
    if (_ifs) {
      _ifs.unget();
      assert(current_state.column > 0);
      --current_state.column;
      current_state.pos = _ifs.tellg();
      current_state.start_symbol_expected = true;
    }
  }
}

size_t DNAFileReader::computeCurrentSequenceLength() {
  if (getCurrentSequenceLength() == size_t(-1)) {
    FileState &current_state = _getState();
    FileState backup_state = current_state;
    gotoSequenceEnd();
    size_t l = current_state.current_sequence_length;
    setState(backup_state);
    current_state.current_sequence_length = l;
  }
  return getCurrentSequenceLength();
}

void DNAFileReader::gotoSequenceEnd(bool check_consistency) {
  switch (getFormat()) {
  case FASTA_FORMAT:
    _gotoSequenceEnd<FASTA_FORMAT>(check_consistency);
    break;
  case FASTQ_FORMAT:
    _gotoSequenceEnd<FASTQ_FORMAT>(check_consistency);
    break;
  default:
    break;
  }
}

void DNAFileReader::gotoSequenceStart() {
  _ifs.clear();
  FileState &current_state = _getState();
  FileReader::setState(current_state.sequence_start_file_state);
  current_state.start_symbol_expected = false;
  current_state.nb_nucleotides = current_state.nb_consecutive_regular_nucleotides = 0;
  current_state.kmer = current_state.kmer_aux = string(current_state.k, '?');
  current_state.kmer_start_pos = 0;
}

bool DNAFileReader::gotoNextSequence(bool check_consistency) {
  gotoSequenceEnd(check_consistency);
  (this->*_parse_mth)();
  return (*this);
}

ostream &operator<<(ostream &os, DNAFileReader::Format fmt) {
  switch (fmt) {
  case DNAFileReader::FASTA_FORMAT:
    os << "Fasta";
    break;
  case DNAFileReader::FASTQ_FORMAT:
    os << "Fastq";
    break;
  default:
    os << "Undefined";
    break;
  }
  return os;
}

bool DNAFileReader::setState(const FileReader::FileState &) {
  Exception e;
  e << "Method " << __FUNCTION__ << "() is not allowed for DNAFileReader instances";
  throw e;
  return false;
}

bool DNAFileReader::setState(const DNAFileReader::FileState &s) {
  FileState &current_state = _getState();
  FileState backup_state = current_state;
  current_state = s;
  bool ok = FileReader::setState(s);
  if (!ok) {
    current_state = backup_state;
  }
  return ok;
}

END_KIM_NAMESPACE
