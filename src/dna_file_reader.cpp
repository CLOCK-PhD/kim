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

DNAFileReader::DNAFileReader(const Settings &settings, const string &filename):
  FileReader(settings, filename)
{}

void DNAFileReader::_onReset() {
  _start_sequence_state = true;
  _nb_nucl = _nb_valid_nucl = 0;
  _current_sequence_description = "";
  _kmer = string(_settings.k(), '?');
}

DNAFileReader::IUPAC DNAFileReader::_toIUPAC(char c) {
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

void DNAFileReader::_processSequenceDescription() {
  _start_sequence_state = false;
  std::getline(_ifs, _current_sequence_description);
  ++_line;
  _col = 0;
  _nb_nucl = _nb_valid_nucl = 0;
}

bool DNAFileReader::_parseSequenceDescription() {
  int c = _nextVisibleCharacter();
  if (c) {
    if (_col == 1) {
      const string &opening_chars = _sequenceStartSymbols();
      if (c != 0) {
        if (opening_chars.find(c) != string::npos) {
          _processSequenceDescription();
        } else {
          ERROR_MSG("Badly formatted file. Character '" << (char) c << "'" << " found while expecting "
                    << ((opening_chars.size() == 1) ? "'" : "one of '") << opening_chars << "'.");
        }
      } else {
        ERROR_MSG("Badly formatted file. End of file found while expecting "
                    << ((opening_chars.size() == 1) ? "'" : "one of '") << opening_chars << "'.");
      }
    } else {
      ERROR_MSG("Badly formatted file. Unexpected character '" << (char) c << "'.");
    }
  }
  return c; // (c != '\0') => Correct Header
}

const string &DNAFileReader::getCurrentKmer() const {
  static const string empty_string;
  return (*this && (_nb_valid_nucl >= _settings.k()) ? _kmer : empty_string);
}

const string &DNAFileReader::getNextKmer() {
  do {
    _parse();
  } while (*this && (_nb_valid_nucl < _settings.k()));
  return getCurrentKmer();
}

const string &DNAFileReader::getForwardKmer(size_t p, bool check_consistency) {
  assert(p >= getCurrentKmerRelativePosition());
  if (check_consistency || (_nb_valid_nucl > p)) {
    while (*this && !_start_sequence_state && (getCurrentKmerRelativePosition() < p)) {
      _parse();
    }
  } else {
    assert(p >= _nb_valid_nucl);
    while (*this && ((_nb_valid_nucl + 1) < p)) {
      _nextVisibleCharacter();
    }
    if (*this) {
      assert((_nb_valid_nucl + 1) == p);
      _parse();
    }
  }
  return getCurrentKmer();
}

bool DNAFileReader::gotoNextSequence(bool check_consistency) {
  bool found = false;
  if (check_consistency) {
    while (*this && !_start_sequence_state) {
      _parse();
    }
    if (*this) {
      found = _parseSequenceDescription();
    }
  } else {
    while (*this && !found) {
      int c = _nextVisibleCharacter();
      if ((_col == 1) and (_sequenceStartSymbols().find(c) != string::npos)) {
        _processSequenceDescription();
        found = true;
      }
    }
  }
  return found;
}

END_KIM_NAMESPACE
