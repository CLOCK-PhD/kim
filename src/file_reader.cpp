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

#include "file_reader.h"

#include "config.h"
#include "kim_exception.h"

#include <cassert>
#include <iostream>
#include <cstdint>

using namespace std;

BEGIN_KIM_NAMESPACE

FileReaderParseError::FileReaderParseError(FileReader &reader):
  Exception()
{
  if (reader.getFilename().empty()) return;
  (*this) << "File '" << reader.getFilename() << "': "
          << "line " << reader.getFileLineNumber() << ", "
          << "column " << reader.getFileColumnNumber()
          << ": ";
}

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


FileReader::FileState::FileState():
  filename(),
  line(-1), column(-1), pos(-1),
  compression_algorithm(FileReader::NO_COMPRESSION),
  current_char(-1), next_char(-1)
{}

FileReader::FileReader(const string &filename, bool warn):
  warn(warn), _state(),
  _ifs_bz2(), _ifs_gzip(), _ifs_none()
{
  assert(_state.filename.empty());
  assert(_state.line == (size_t) -1);
  assert(_state.column == (size_t) -1);
  assert(_state.pos == -1);
  assert(_state.compression_algorithm == CompressionAlgorithm::NO_COMPRESSION);
  assert(_state.current_char == -1);
  assert(_state.next_char == -1);
  assert(!_ifs_bz2.is_open());
  assert(!_ifs_gzip.is_open());
  assert(!_ifs_none.is_open());
  if (!filename.empty()){
    open(filename);
  }
}

FileReader::~FileReader() {
  close();
}

bool FileReader::good() const {
  bool ok = false;
  switch (getState().compression_algorithm) {
  case CompressionAlgorithm::NO_COMPRESSION:
    assert(!ok);
    assert(!_ifs_gzip.is_open());
    assert(!_ifs_bz2.is_open());
    ok = _ifs_none.good();
    break;
  case CompressionAlgorithm::GZIP:
    assert(!ok);
    assert(!_ifs_none.is_open());
    assert(!_ifs_bz2.is_open());
    ok = _ifs_gzip.good();
    break;
  case CompressionAlgorithm::BZ2:
    assert(!ok);
    assert(!_ifs_none.is_open());
    assert(!_ifs_gzip.is_open());
    ok = _ifs_bz2.good();
    break;
  }
  return ok;
}

void FileReader::clear() {
  FileState &current_state = _getState();
  switch (current_state.compression_algorithm) {
  case CompressionAlgorithm::NO_COMPRESSION:
    _ifs_none.clear();
    break;
  case CompressionAlgorithm::GZIP:
    _ifs_gzip.clear();
    break;
  case CompressionAlgorithm::BZ2:
    _ifs_bz2.clear();
    break;
  }
}

void FileReader::open(const string &filename) {
  close();
  assert(!_ifs_bz2.is_open());
  assert(!_ifs_gzip.is_open());
  assert(!_ifs_none.is_open());
  FileState &current_state = _getState();
  current_state.filename = filename;
  _ifs_none.open(current_state.filename.c_str());
  if (*this) {
    current_state.line = current_state.column = 0;
    _detectCompressionAlgorithm();
    _onOpen();
  } else {
    current_state.filename = "";
    if (warn) {
      cerr << "Unable to open file '" << filename << "'." << endl;
    }
  }
}

void FileReader::close() {
  if (_ifs_none.is_open()) _ifs_none.close();
  if (_ifs_gzip.is_open()) _ifs_gzip.close();
  if (_ifs_bz2.is_open()) _ifs_bz2.close();
  assert(!_ifs_none.is_open());
  assert(!_ifs_gzip.is_open());
  assert(!_ifs_bz2.is_open());
  FileState &current_state = _getState();
  current_state = FileState();
  reset();
  _onClose();
  assert(!_ifs_none.is_open());
  assert(!_ifs_gzip.is_open());
  assert(!_ifs_bz2.is_open());
}

void FileReader::reset() {
  FileState &current_state = _getState();
  bool is_open = false;
  clear();
  _seekg(0);
  is_open = _is_open();
  current_state.current_char = -1;
  current_state.next_char = -1;
  current_state.pos = (is_open ? 0 : -1);
  current_state.line = current_state.column = (is_open ? 0 : -1);
  _onReset();
}

bool FileReader::setState(const FileState &s) {
  bool ok = true;
  FileState &current_state = _getState();
  if (s.filename != current_state.filename) {
    close();
    if (!s.filename.empty()) {
      open(s.filename);
    }
  }

  ok &= (s.filename == current_state.filename);
  ok &= (s.compression_algorithm == current_state.compression_algorithm);
  if (ok) {
    current_state.line = s.line;
    current_state.column = s.column;
    clear();
#ifndef NDEBUG
    ifstream::pos_type p =
#endif
      _seekg(s.pos);
    assert(p == s.pos);
    current_state.pos = s.pos;
    current_state.current_char = s.current_char;
    current_state.next_char = s.next_char;
    if (current_state.next_char != -1) {
      _loadNextChar();
    }
  }
  return ok;
}

void FileReader::_detectCompressionAlgorithm() {
  if (!_ifs_none.is_open()) return;
  FileState &current_state = _getState();
  assert(current_state.line == 0);
  assert(current_state.column == 0);
  unsigned char magic[2] = { 0, 0 };
  _ifs_none.read((char *)magic, sizeof(magic));
  current_state.compression_algorithm = CompressionAlgorithm::NO_COMPRESSION;
  if ((magic[0] == 0x1f) && (magic[1] == 0x8b)) {
    current_state.compression_algorithm = CompressionAlgorithm::GZIP;
    assert(!_ifs_gzip.is_open());
    _ifs_gzip.open(current_state.filename.c_str());
    assert(_ifs_gzip.is_open());
    assert(_ifs_gzip.good());
    _ifs_none.close();
    assert(!_ifs_none.is_open());
  } else if ((magic[0] == 0x42) && (magic[1] == 0x5a)) {
    current_state.compression_algorithm = CompressionAlgorithm::BZ2;
    assert(!_ifs_bz2.is_open());
    _ifs_bz2.open(current_state.filename.c_str());
    assert(_ifs_bz2.is_open());
    assert(_ifs_bz2.good());
    _ifs_none.close();
    assert(!_ifs_none.is_open());
  }
  reset();
}

void FileReader::_loadNextChar() {
  assert(good());
  FileState &current_state = _getState();
  switch (current_state.compression_algorithm) {
  case CompressionAlgorithm::NO_COMPRESSION:
    assert(_ifs_none.is_open());
    current_state.next_char = _ifs_none.get();
   break;
  case CompressionAlgorithm::BZ2:
    assert(_ifs_bz2.is_open());
    current_state.next_char = _ifs_bz2.get();
    break;
  case CompressionAlgorithm::GZIP:
    assert(_ifs_gzip.is_open());
    current_state.next_char = _ifs_gzip.get();
    break;
  }
}

int FileReader::peek() {
  FileState &current_state = _getState();
  if ((current_state.next_char == -1) && good()) {
    _loadNextChar();
  }
  return current_state.next_char;
}

int FileReader::get() {
  FileState &current_state = _getState();
  switch (current_state.current_char) {
  case '\n':
    ++current_state.line;
    current_state.column = 0;
    /* FALLTHROUGH */
  case -1:
    break;
  default:
    ++current_state.column;
  }
  current_state.current_char = peek();
  current_state.next_char = -1;
  current_state.pos += 1;

  return current_state.current_char;
}


string FileReader::getline(const char delim) {
  string res;
  int c;
  while (good() && ((c = get()) != delim)) {
    if (c != -1) {
      res += c;
    }
  }
  return res;
}

size_t FileReader::ignore(const char delim) {
  FileState &current_state = _getState();
  ifstream::pos_type p = current_state.pos;
  while (good() && (get() != delim));
  return current_state.pos - p + 1;
}

ifstream::pos_type FileReader::_seekg(ifstream::pos_type p) {
  FileState &current_state = _getState();
  switch (current_state.compression_algorithm) {
  case CompressionAlgorithm::NO_COMPRESSION:
    _ifs_none.seekg(p);
    break;
  case CompressionAlgorithm::BZ2:
    _ifs_bz2.seekg(p);
    break;
  case CompressionAlgorithm::GZIP:
    _ifs_gzip.seekg(p);
    break;
  }
  return p;
}

ifstream::pos_type FileReader::_tellg() {
  FileState &current_state = _getState();
  return current_state.pos;
}

char FileReader::_nextVisibleCharacter(bool stop_before) {
  int c = -1;
  while (*this && (peek() <= 32)) {
    c = get();
    switch (c) {
    case '\t':
    case ' ':
    case '\n':
    case -1:
      break;
    default:
      WARNING_MSG("Unexpected non printable character (code " << c << ")");
    }
  }

  if (stop_before) {
    c = peek();
  } else {
    c = get();
  }
  return ((c == -1) ? 0 : c);
}

const char *FileReader::_search_directories[] = {
                                                 "./",
                                                 PACKAGE_DATADIR "/",
                                                 /* The following directories are mostly for development purpose */
                                                 "resources/",
                                                 "../resources/",
                                                 "../",
                                                 /* The last array entry must be NULL */
                                                 NULL
};

string FileReader::findFile(const string &filename, const char **directories) {
  if (filename.empty()) return filename;
  string res = filename;
  const char *dir = NULL;
  do {
    ifstream ifs(res);
    if (ifs) {
      ifs.close();
      dir = NULL;
    } else {
      if (directories && *directories && (filename[0] != '/')) {
        dir = *directories++;
        res = dir;
        if (res.back() != '/') {
          res += '/';
        }
        res += filename;
      } else {
        dir = NULL;
        res = "";
      }
    }
  } while (dir);
  return res;
}

END_KIM_NAMESPACE
