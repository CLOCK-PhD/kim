/******************************************************************************
*                                                                             *
*  Copyright © 2023      -- IGH / LIRMM / CNRS / UM                           *
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

#include "file_reader.h"
#include "config.h"

#include <sstream>
#include <iostream>
#include <exception>

using namespace std;

BEGIN_KIM_NAMESPACE


/////////////// Some stuff to handling parse error ///////////////

class FileReaderParseError: public exception {

private:
  string s;

public:

  FileReaderParseError(FileReader &reader) {
    if (reader.getFilename().empty()) return;
    stringstream ss;
    ss << "File '" << reader.getFilename() << "'"
       << ": line " << reader.getFileLineNumber()
       << ": column " << reader.getFileColumnNumber()
       << ": ";
    s = ss.str();
  }

  inline virtual const char *what() const noexcept {
    return s.c_str();
  }

  template <typename T>
  inline FileReaderParseError &operator<<(const T &t) {
    stringstream ss;
    ss << t;
    s += ss.str();
    return *this ;
  }

};

#define ALERT_MSG(header, msg)  \
  if (warn) {                   \
    cerr << header << ":"       \
         << filename << ":"     \
         << (line + 1) << ":"   \
         << col << ": "         \
         << msg << endl;        \
  }                             \
  (void) 0

#define WARNING_MSG(msg) ALERT_MSG("Warning", msg)
#define ERROR_MSG(msg) FileReaderParseError e(*this); e << msg; throw e

//////////////////////////////////////////////////////////////////


FileReader::FileReader(const char *filename, size_t k, bool warn) {
  open(filename, k, warn);
}

FileReader::~FileReader() {
  ifs.close();
}

void FileReader::open(const char *filename, size_t k, bool warn) {
  if (ifs.is_open()) {
    ifs.close();
  }
  if (!k) {
    ERROR_MSG("The length of k-mer must be strictly positive!");
  }
  this->filename = filename;
  ifs.open(filename);
  if (ifs) {
    line = col = 0;
    start_sequence_state = true;
    nb_nucl = nb_valid_nucl = 0;
    this->k = k;
    read_id = "";
    kmer = string(k, '?');
    kmer_pos = 0;
    this->warn = warn;
  }
}

const string &FileReader::getFilename() const {
  return filename;
}

size_t FileReader::getFileLineNumber() const {
  return line + 1;
}

size_t FileReader::getFileColumnNumber() const {
  return col + 1;
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

char FileReader::nextVisibleCharacter() {
  int c = -1;
  while (ifs && ((c = ifs.get()) <= 32)) {
    switch (c) {
    case '\t':
    case ' ':
      ++col;
      break;
    case '\n':
      ++line;
      col = 0;
      break;
    default:
      if (ifs) {
        WARNING_MSG("Unexpected non printable character (code " << c << ")");
      }
    }
  }
  ++col;
  return ((c != -1) ? c : 0);
}

const string &FileReader::getNextKmer() {
  do {
    int c = nextVisibleCharacter();
    if (c) {
      if (start_sequence_state) {
        if (col == 1) {
          switch (c) {
          case '@':
            start_sequence_state = false;
            std::getline(ifs, read_id);
            ++line;
            col = 0;
            nb_nucl = 0;
            nb_valid_nucl = 0;
            kmer_pos = 0;
            break;
          default:
            ERROR_MSG("Badly formatted file. Character '" << (char) c << "' found while expecting '@'.");
          }
        } else {
          ERROR_MSG("Badly formatted file. Unexpected character '" << (char) c << "'.");
        }
      } else {
        uint8_t nucl = getIUPACNucleotide(c);
        switch (nucl) {
        case IUPAC_A:
        case IUPAC_C:
        case IUPAC_G:
        case IUPAC_T:
          if (nb_valid_nucl >= k) {
            for (size_t i = 1; i < k; ++i) {
              kmer[i - 1] = kmer[i];
            }
            kmer[k - 1] = c;
          } else {
            kmer[nb_valid_nucl] = c;
          }
          if (++nb_valid_nucl >= k) {
            ++kmer_pos;
          }
          ++nb_nucl;
          break;
        case IUPAC_GAP:
          // ignore gaps
          break;
        case IUPAC_UNDEFINED:
          if ((col == 1) && (c == '+')) {
            // End of the sequence, let's skip the quality
            string header;
            std::getline(ifs, header);
            ++line;
            col = 0;
            if (!header.empty() && (header != read_id)) {
              WARNING_MSG("Badly formatted file. The header following the '+' sign (" << header << ") is not the expected one (" << read_id << ").");
            }
            while (nb_nucl--) {
              c = nextVisibleCharacter();
            }
            start_sequence_state = true;
          } else {
            WARNING_MSG("Badly formatted file. Unexpected character '" << (char) c << "'.");
          }
        default:
          // This is a degeneracy symbol
          ++nb_nucl;
          nb_valid_nucl = 0;
        }
      }
    } else {
      kmer_pos = 0;
    }
  } while (ifs && !kmer_pos);
  return getCurrentKmer();
}

const string &FileReader::getCurrentKmer() const {
  static const string empty_string;
  return (kmer_pos ? kmer : empty_string);
}

const string &FileReader::getCurrentRead() const {
  return read_id;
}

size_t FileReader::getCurrentKmerRelativePosition() const {
  return kmer_pos;
}

END_KIM_NAMESPACE
