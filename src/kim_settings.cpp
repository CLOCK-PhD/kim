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

#include "kim_settings.h"

#include "config.h"

#include <iostream>

using namespace std;

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

Settings::Settings(size_t k, size_t p, bool warn, bool freeze):
  _k(0), _p(0), _s(0),
  _warn(false), _frozen(false) {
  if (k || p) {
    setKmerLength(k);
    _warn = warn;
    setKmerPrefixLength(p);
  }
  _warn = warn;
  _frozen = freeze;
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
  if (_p + 1 >= _k) {
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

void Settings::warn(bool status) {
  CHECK_FROZEN_STATE(!frozen(), warn);
  _warn = status;
}

END_KIM_NAMESPACE
