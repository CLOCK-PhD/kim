/******************************************************************************
*                                                                             *
*  Copyright © 2025      -- IGH / LIRMM / CNRS / UM                           *
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

#include "variant_allele_frequencies.h"

#include "config.h"

#include <cassert>

using namespace std;

BEGIN_KIM_NAMESPACE

VariantAlleleFrequencies::VariantPopulationAlleleFrequency::VariantPopulationAlleleFrequency(const string &v, const map<string, float> &population_af):
  variant(v), population_af(population_af) {
}

VariantAlleleFrequencies::VariantAlleleFrequencies(const set<string> &tags):
  tags(tags) {
}

VariantAlleleFrequencies &VariantAlleleFrequencies::compute(const VariantKmerEnumerator &vke) {

  _vpaf_array.clear();

  for (auto const &v: vke.getCurrentVariantIDs()) {
    _vpaf_array.emplace_back(v);
  }

  for (auto const &tag: tags) {
    vector<float> v_af;
    if (vke.getAlleleFrequencies(v_af, tag)) {
      size_t i = 0;
      for (auto &vpaf: _vpaf_array) {
        vpaf.population_af[tag] = v_af[i++];
      }
    }
  }

  if (_vpaf_array.empty() || _vpaf_array[0].population_af.empty()) {
    _vpaf_array.clear();
  }

  return *this;

}

vector<string> VariantAlleleFrequencies::getAlternateAlleleVariants() const {
  vector<string> names;
  names.reserve(_vpaf_array.size());
  for (auto const &vpaf: _vpaf_array) {
    names.push_back(vpaf.variant);
  }
  return names;
}

const VariantAlleleFrequencies::VariantPopulationAlleleFrequency &VariantAlleleFrequencies::getAlleleFrequency(const string &alternate_allele) const {
  size_t a = 0, b = _vpaf_array.size(), m = 0;
  bool ok = false;
  while (!ok && (a < b)) {
    m = (a + b) / 2;
    const int v = _vpaf_array[m].variant.compare(alternate_allele);
    if (v < 0) {
      a = m + 1;
    } else if (v > 0) {
      b = m;
    } else {
      ok = true;
    }
  }
  if (!ok) {
    Exception e;
    e << "There is no variant with ID '" << alternate_allele << "'";
    throw e;
  }
  return _vpaf_array[m];
}


void VariantAlleleFrequencies::printCSVHeader(ostream &os) const {
  os << "#Variant";
  for (auto const &tag: tags) {
    os << "\t" << tag;
  }
  os << "\n";
}

void VariantAlleleFrequencies::toCSV(ostream &os) const {
  for (auto const &vpaf: _vpaf_array) {
    os << vpaf.variant;
    for (auto const &tag: tags) {
      map<string, float>::const_iterator it = vpaf.population_af.find(tag);
      os << "\t" << (it == vpaf.population_af.cend() ? float(-1) : it->second);
    }
    os << "\n";
  }
}

ostream &operator<<(ostream& os, const VariantAlleleFrequencies &vaf) {
  vaf.toCSV(os);
  return os;
}

END_KIM_NAMESPACE
