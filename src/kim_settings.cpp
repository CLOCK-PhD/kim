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

#include "kim_settings.h"

#include "config.h"

#include <cassert>
#include <iostream>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

BEGIN_KIM_NAMESPACE

#define ERROR_MSG(msg)          \
  do {                          \
    BadSettingsException error; \
    error << msg;               \
    throw error;                \
  } while (0)

#define CHECK_FROZEN_STATE(expected_state, mth)                         \
  if (!(expected_state)) {                                              \
    ERROR_MSG("Settings must be " << (expected_state ? "frozen" : "unfrozen") \
              << " before calling Settings::" << #mth                   \
              << "() method.");                                         \
  }                                                                     \
  (void) 0

Settings::Settings(size_t k, size_t p, const fs::path &index_directory,
                   bool warn, bool check_consistency, bool allow_overwrite,
                   double alpha, double threshold, bool weak_mode,
                   bool freeze):
  _k(0), _p(0), _s(0), _index_directory(index_directory),
  _warn(false), _check_consistency(check_consistency),
  _allow_overwrite(allow_overwrite),
  _alpha(alpha), _threshold(threshold), _weak_mode(weak_mode),
  _frozen(false) {
  if (k || p) {
    setKmerLength(k);
    _warn = warn;
    setKmerPrefixLength(p);
  }
  _warn = warn;
  _frozen = freeze;
  assert(_alpha >= 0);
  assert(_alpha <= 1);
  assert(_threshold >= 0);
  assert(_threshold <= 1);
}

void Settings::freeze() {
  if (!valid()) {
    ERROR_MSG("Settings aren't valid, thus calling the Settings::freeze() method is not allowed.");
  }
  _frozen = true;
}

void Settings::setKmerLength(size_t k) {
  CHECK_FROZEN_STATE(!frozen(), setKmerLength);
  if (k < 2) {
    ERROR_MSG("Can't set the length of the k-mers to " << k
              << " (length must be greater or equal to 2).");
  }
  _k = k;
  if (_p + 1 >= _k) {
    if (_warn) {
      cerr << "The length of k-mers is set to " << _k
           << " but current length of the prefixes was " << _p << "."
           << " Setting the length of the prefixes to " << (_k - 1)
           << "." << endl;
    }
    _p = _k - 1;
  }
  _s = _k - _p;
}

void Settings::setKmerPrefixLength(size_t p) {
  CHECK_FROZEN_STATE(!frozen(), setKmerPrefixLength);
  if (p == 0) {
    ERROR_MSG("Can't set the prefix length of the k-mers to 0 (it must be a strictly positive value).");
  }
  _p = p;
  if (_p + 1 > _k) {
    if (_warn) {
      cerr << "The length of k-mers is set to " << _k
           << " but wanted length of the prefixes is " << _p << "."
           << " Setting the length of the prefixes to " << (_k - 1)
           << "." << endl;
    }
    _p = _k - 1;
  }
  _s = _k - _p;
}

void Settings::validateDirectory(const fs::path &path, bool must_exist, bool must_not_exist) {
  if (must_exist || must_not_exist) {
    fs::file_status s = fs::status(path);
    if (must_exist) {
      if (!fs::is_directory(s)) {
        BadSettingsException e;
        e << "The directory '" << path << "' doesn't exist or you don't have enough access rights";
        throw e;
      }

      fs::perms permissions = s.permissions();
      static const fs::perms urx = fs::perms::owner_read | fs::perms::owner_exec;
      static const fs::perms grx = fs::perms::group_read | fs::perms::group_exec;
      static const fs::perms orx = fs::perms::others_read | fs::perms::others_exec;
      if (((permissions & urx) != urx)       // permissions doesn't match 'dr.x......'
          && ((permissions & grx) != grx)    // permissions doesn't match 'd...r.x...'
          && ((permissions & orx) != orx)) { // permissions doesn't match 'd......r.x'
        BadSettingsException e;
        e << "The directory '" << path << "' is not readable";
        throw e;
      }
    }

    if (must_not_exist && fs::exists(s)) {
      BadSettingsException e;
      e << "The directory '" << path << "' already exists";
      throw e;
    }
  }
}

void Settings::setIndexDirectory(const fs::path &path, bool must_exist, bool must_not_exist) {
  CHECK_FROZEN_STATE(!frozen(), setIndexDirectory);
  validateDirectory(path, must_exist, must_not_exist);
  _index_directory = path;
}

void Settings::warn(bool status) {
  CHECK_FROZEN_STATE(!frozen(), warn);
  _warn = status;
}

void Settings::checkConsistency(bool status) {
  CHECK_FROZEN_STATE(!frozen(), warn);
  _check_consistency = status;
}

void Settings::allowOverwrite(bool status) {
  CHECK_FROZEN_STATE(!frozen(), warn);
  _allow_overwrite = status;
}

void Settings::alpha(double v) {
  CHECK_FROZEN_STATE(!frozen(), alpha);
  assert(v >= 0);
  assert(v <= 1);
  _alpha = v;
}

void Settings::threshold(double v) {
  CHECK_FROZEN_STATE(!frozen(), threshold);
  assert(v >= 0);
  assert(v <= 1);
  _threshold = v;
}

void Settings::weakMode(bool b) {
  CHECK_FROZEN_STATE(!frozen(), weakMode);
  _weak_mode = b;
  assert(weakMode() == !strictMode());
}

void Settings::strictMode(bool b) {
  CHECK_FROZEN_STATE(!frozen(), strictMode);
  _weak_mode = !b;
  assert(strictMode() == !weakMode());
}

bool Settings::addAlleleFrequencyTag(const string &tag) {
  CHECK_FROZEN_STATE(!frozen(), addAlleleFrequencyTag);
  return _allele_frequency_tags.insert(tag).second;
}

bool Settings::removeAlleleFrequencyTag(const string &tag) {
  CHECK_FROZEN_STATE(!frozen(), removeAlleleFrequencyTag);
  return _allele_frequency_tags.erase(tag) == 1;
}

END_KIM_NAMESPACE
