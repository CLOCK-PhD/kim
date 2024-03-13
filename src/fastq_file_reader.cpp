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

#include "fastq_file_reader.h"

#include "config.h"
#include "kim_exception.h"

#include <cassert>
#include <iostream>

using namespace std;

BEGIN_KIM_NAMESPACE

#define WARNING_MSG(msg)          \
  if (_settings.warn()) {         \
    cerr << "Warning:"            \
         << _filename << ":"      \
         << (_line + 1) << ":"    \
         << _col << ": "          \
         << msg << endl;          \
  }                               \
  (void) 0

#define ERROR_MSG(msg)                          \
  do {                                          \
    FileReaderParseError error(*this);          \
    error << msg;                               \
    throw error;                                \
  } while (0)

FastqFileReader::FastqFileReader(const Settings &settings): FileReader(settings) {}

FastqFileReader::FastqFileReader(const char *filename, const Settings &settings):
  FileReader(settings)
{
  open(filename);
}

FastqFileReader::FastqFileReader(const string &filename, const Settings &settings):
  FileReader(settings)
{
  open(filename);
}

void FastqFileReader::_onReset() {
  _start_sequence_state = true;
  _nb_nucl = _nb_valid_nucl = 0;
  _read_id = "";
  _kmer = string(_settings.k(), '?');
}

#define IUPAC_A 1
#define IUPAC_C 2
#define IUPAC_G 4
#define IUPAC_T 8
#define IUPAC_GAP 16
#define IUPAC_UNDEFINED 0

uint8_t getIUPACNucleotide(char c) {
  if ((c >= 'a') && (c <= 'z')) {
    c += 'A' - 'a'; // Upcase
  }
  switch (c) {
  case 'A': return IUPAC_A;
  case 'C': return IUPAC_C;
  case 'G': return IUPAC_G;
  case 'T':
  case 'U': return IUPAC_T;
  case 'R': return IUPAC_A | IUPAC_G; // Purine
  case 'Y': return IUPAC_C | IUPAC_T; // Pyrimidine
  case 'K': return IUPAC_G | IUPAC_T; // Ketones
  case 'M': return IUPAC_A | IUPAC_C; // Amino groups
  case 'S': return IUPAC_C | IUPAC_G; // Strong
  case 'W': return IUPAC_A | IUPAC_T; // Weak
  case 'B': return IUPAC_C | IUPAC_G | IUPAC_T; // Not A
  case 'D': return IUPAC_A | IUPAC_G | IUPAC_T; // Not C
  case 'H': return IUPAC_A | IUPAC_C | IUPAC_T; // Not G
  case 'V': return IUPAC_A | IUPAC_C | IUPAC_G; // Not T neither U
  case 'N': return IUPAC_A | IUPAC_C | IUPAC_G | IUPAC_T; // Any
  case '.':
  case '-': return IUPAC_GAP; // Gap is valid in IUPAC even if it's
                              // not a nucleotide
  default:
    return IUPAC_UNDEFINED;
  }
}

void FastqFileReader::_processSequenceHeader() {
  _start_sequence_state = false;
  std::getline(_ifs, _read_id);
  ++_line;
  _col = 0;
  _nb_nucl = _nb_valid_nucl = 0;
}

bool FastqFileReader::_parseSequenceHeader() {
  int c = _nextVisibleCharacter();
  if (c) {
    if (_col == 1) {
      switch (c) {
      case '@':
        _processSequenceHeader();
        break;
      case 0:
        ERROR_MSG("Badly formatted file. End of file found while expecting '@'.");
      default:
        ERROR_MSG("Badly formatted file. Character '" << (char) c << "' found while expecting '@'.");
      }
    } else {
      ERROR_MSG("Badly formatted file. Unexpected character '" << (char) c << "'.");
    }
  }
  return c;
}

void FastqFileReader::_parseQualitySeparator() {
  // End of the sequence, let's skip the quality
  string header;
  std::getline(_ifs, header);
  ++_line;
  _col = 0;
  if (!header.empty() && (header != _read_id)) {
    WARNING_MSG("Badly formatted file. The header following the '+' sign (" << header << ") is not the expected one (" << _read_id << ").");
  }
}

void FastqFileReader::_parseQualitySequence() {
  int c = -1;
  while (c && _nb_nucl--) {
    c = _nextVisibleCharacter();
  }
  if (!c) {
    ERROR_MSG("Badly formatted file. End of file found while " << _nb_nucl << " quality symbols were still expected.");
  }
  _start_sequence_state = true;
}

void FastqFileReader::_parse() {
  bool ok;
  do {
    ok = true;
    if (_start_sequence_state) {
      ok = !_parseSequenceHeader();
    } else {
      int c = _nextVisibleCharacter();
      if (c) {
        uint8_t nucl = getIUPACNucleotide(c);
        switch (nucl) {
        case IUPAC_A:
        case IUPAC_C:
        case IUPAC_G:
        case IUPAC_T:
          if (_nb_valid_nucl >= _settings.k()) {
            for (size_t i = 1; i < _settings.k(); ++i) {
              _kmer[i - 1] = _kmer[i];
            }
            _kmer[_settings.k() - 1] = c;
          } else {
            _kmer[_nb_valid_nucl] = c;
          }
          ++_nb_nucl;
          ++_nb_valid_nucl;
          break;
        case IUPAC_GAP:
          // ignore gaps
          ok = false;
          break;
        case IUPAC_UNDEFINED:
          if ((_col == 1) && (c == '+')) {
            _parseQualitySeparator();
            _parseQualitySequence();
          } else {
            WARNING_MSG("Badly formatted file. Unexpected character '" << (char) c << "'.");
          }
          ok = false;
          break;
        default:
          // This is a degeneracy symbol
          ++_nb_nucl;
          _nb_valid_nucl = 0;
        }
      } else {
        _nb_nucl = _nb_valid_nucl = 0;
      }
    }
  } while (*this && !ok);
}

const string &FastqFileReader::getNextKmer() {
  do {
    _parse();
  } while (*this && (_nb_valid_nucl < _settings.k()));
  return getCurrentKmer();
}

const string &FastqFileReader::getCurrentKmer() const {
  static const string empty_string;
  return (*this && (_nb_valid_nucl >= _settings.k()) ? _kmer : empty_string);
}

bool FastqFileReader::gotoNextSequence(bool check_consistency) {
  bool found = false;
  if (check_consistency) {
    while (*this && !_start_sequence_state) {
      int c = _nextVisibleCharacter();
      if (c) {
        uint8_t nucl = getIUPACNucleotide(c);
        switch (nucl) {
        case IUPAC_GAP:
          // cout << "=> A gap" << endl;
          // ignore gaps
          break;
        case IUPAC_UNDEFINED:
          if ((_col == 1) && (c == '+')) {
            _parseQualitySeparator();
            _parseQualitySequence();
          } else {
            WARNING_MSG("Badly formatted file. Unexpected character '" << (char) c << "'.");
          }
          break;
        default:
          // This is a nucleotide symbol (degeneracy or not)
          ++_nb_nucl;
        }
      }
    }
    if (*this) {
      found = _parseSequenceHeader();
    }
  } else {
    while (*this && !found) {
      int c = _nextVisibleCharacter();
      if ((_col == 1) and (c == '@')) {
        _processSequenceHeader();
        found = true;
      } else {
        ++_col;
      }
    }
  }
  return found;
}

END_KIM_NAMESPACE
