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

FastqFileReader::FastqFileReader(const Settings &settings, const string &filename):
  DNAFileReader(settings, filename)
{}

static const string _FASTQ_START_SEQUENCE_SYMBOLS = "@";
const string &FastqFileReader::_sequenceStartSymbols() const {
  return _FASTQ_START_SEQUENCE_SYMBOLS;
}


void FastqFileReader::_parseQualitySeparator() {
  // End of the sequence, let's skip the quality
  string description;
  std::getline(_ifs, description);
  ++_line;
  _col = 0;
  if (!description.empty() && (description != _current_sequence_description)) {
    WARNING_MSG("Badly formatted file. The description following the '+' sign (" << description << ") is not the expected one (" << _current_sequence_description << ").");
  }
}

void FastqFileReader::_parseQualitySequence() {
  int c = -1;
  size_t nb = _nb_nucl;
  while (c && nb--) {
    c = _nextVisibleCharacter();
  }
  if (!c) {
    ERROR_MSG("Badly formatted file. End of file found while " << nb << " quality symbols were still expected.");
  }
  _start_sequence_state = true;
}

void FastqFileReader::_parse() {
  bool ok;
  if (_start_sequence_state) {
    ok = _parseSequenceDescription();
  } else {
    do {
      int c = _nextVisibleCharacter();
      if (c) {
        IUPAC nucl = _toIUPAC(c);
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
          ok = (_nb_valid_nucl >= _settings.k());
          break;
        case IUPAC_GAP:
          // ignore gaps
          ok = false;
          break;
        case IUPAC_UNDEFINED:
          if ((_col == 1) && (c == '+')) {
            _parseQualitySeparator();
            _parseQualitySequence();
            _nb_valid_nucl = 0;
            ok = true;
          } else {
            WARNING_MSG("Badly formatted file. Unexpected character '" << (char) c << "'.");
          }
          break;
        default:
          // This is a degeneracy symbol
          ++_nb_nucl;
          _nb_valid_nucl = 0;
        }
      } else {
        assert(!*this);
        ok = false;
        _nb_nucl = _nb_valid_nucl = 0;
      }
    } while (*this && !ok);
  }
  assert(ok xor !*this);
}

END_KIM_NAMESPACE
