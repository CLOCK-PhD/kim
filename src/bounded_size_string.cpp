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

#include "bounded_size_string.h"

#include "config.h"

#include <cstring>

using namespace std;

BEGIN_KIM_NAMESPACE

size_t BoundedSizeString::_maximal_size = 0;

size_t BoundedSizeString::_nb_instances = 0;

bool BoundedSizeString::setMaximalSize(size_t maximal_size) {
  bool res = (_nb_instances == 0);
  if (res) {
    _maximal_size = maximal_size;
  }
  return res;
}

void BoundedSizeString::_copy(const char *s) {
  if (!_str) {
    _str = new char[_maximal_size + 1];
    ++_nb_instances;
  }
  size_t n = s ? strlen(s) : 0;
  if (n > _maximal_size) {
    n = _maximal_size;
  }
  memcpy(_str, s, n);
  memset(_str + n, 0, _maximal_size + 1 - n);
}

BoundedSizeString::BoundedSizeString(const char *c_str):
  _str(NULL) {
  _copy(c_str);
}

BoundedSizeString::BoundedSizeString(const string &s):
  _str(NULL) {
  _copy(s.c_str());
}

BoundedSizeString::BoundedSizeString(const BoundedSizeString &s):
  _str(NULL) {
  _copy(s._str);
}

BoundedSizeString::~BoundedSizeString() {
  --_nb_instances;
  delete [] _str;
}

BoundedSizeString &BoundedSizeString::operator=(const BoundedSizeString &s) {
  if (this != &s) {
    _copy(s._str);
  }
  return *this;
}

BoundedSizeString &BoundedSizeString::operator=(const string &s) {
  if (_str != s.c_str()) {
    _copy(s.c_str());
  }
  return *this;
}

BoundedSizeString &BoundedSizeString::operator=(const char *s) {
  if (_str != s) {
    _copy(s);
  }
  return *this;
}

int BoundedSizeString::compare(const BoundedSizeString &s) const {
  if (!_maximal_size) return 0;
  size_t i = 0;
  while ((_str[i] != '\0') && (s._str[i] != '\0') && (_str[i] == s._str[i])) {
    ++i;
  }
  return ((_str[i] > s._str[i])
          ? 1
          : ((_str[i] == s._str[i])
             ? 0
             : -1));
}

int BoundedSizeString::reverse_compare(const BoundedSizeString &s) const {
  if (!_maximal_size) return 0;
  size_t n1 = length();
  size_t n2 = s.length();
  if (n1 < n2) return -1;
  if (n1 > n2) return 1;
  size_t i = n1;
  do {
    --i;
  } while (i && (_str[i] == s._str[i]));
  return ((_str[i] > s._str[i])
          ? 1
          : ((_str[i] == s._str[i])
             ? 0
             : -1));
}

size_t BoundedSizeString::length() const {
  return strlen(_str);
}

void BoundedSizeString::clear() {
  // Long and proper version:
  memset(_str, 0, _maximal_size + 1);
  // Faster and dirty version:
  // if (_str) {
  //   _str[0] = '\0';
  // }
}

bool BoundedSizeString::empty() const {
  return (!_str || _str[0] == '\0');
}

istream &getline(istream &is, BoundedSizeString& s, char delim) {
  size_t i = 0;
  int c;
  while (is && (i < BoundedSizeString::getMaximalSize()) && ((c = is.get()) != delim)) {
    s[i++] = c;
  }
  memset(&s[i], 0, BoundedSizeString::getMaximalSize() + 1 - i);
  return is;
}

istream &operator>>(istream &is, BoundedSizeString &s) {
  size_t i = 0;
  while (is && (is.peek() <= ' ')) {
    is.get();
  }
  while (is && (i < BoundedSizeString::getMaximalSize()) && (is.peek() > ' ')) {
    s[i++] = is.get();
  }
  memset(&s[i], 0, BoundedSizeString::getMaximalSize() + 1 - i);
  return is;
}

END_KIM_NAMESPACE
