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

#include "file_reader.h"

#include "config.h"
#include "kim_exception.h"

#include <cassert>
#include <iostream>

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

FileReader::FileReader(const string &filename, bool warn): warn(warn) {
  if (!filename.empty()){
    open(filename);
  }
}

FileReader::~FileReader() {
  close();
}

void FileReader::open(const string &filename) {
  close();
  FileState &current_state = _getState();
  current_state.filename = filename;
  _ifs.open(current_state.filename.c_str());
  if (!*this) {
    current_state.filename = "";
    if (warn) {
      cerr << "Unable to open file '" << filename << "'." << endl;
    }
  }
  _onOpen();
}

void FileReader::close() {
  if (_ifs.is_open()) {
    _ifs.close();
  }
  _getState().filename = "";
  reset();
  _onClose();
}

void FileReader::reset() {
  _ifs.clear();
  if (_ifs) {
    _ifs.seekg(0);
  }
  FileState &current_state = _getState();
  current_state.pos = _ifs.tellg();
  current_state.line = current_state.column = 0;
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
  ok = (s.filename == current_state.filename);
  if (ok) {
    _ifs.clear();
    current_state.line = s.line;
    current_state.column = s.column;
    _ifs.seekg(s.pos);
    current_state.pos = _ifs.tellg();
  }
  return ok;
}

char FileReader::_nextVisibleCharacter() {
  int c = -1;
  FileState &current_state = _getState();
  while (*this && ((c = _ifs.get()) <= 32)) {
    switch (c) {
    case '\t':
    case ' ':
      ++current_state.column;
      break;
    case '\n':
      ++current_state.line;
      current_state.column = 0;
      break;
    default:
      if (*this) {
        WARNING_MSG("Unexpected non printable character (code " << c << ")");
      }
    }
  }
  ++current_state.column;
  current_state.pos = _ifs.tellg();
  return ((c != -1) ? c : 0);
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
