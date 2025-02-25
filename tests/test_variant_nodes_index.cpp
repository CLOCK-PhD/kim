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

#include <variant_nodes_index.h>

#include <string>
#include <iostream>
#ifdef NDEBUG
#  undef NDEBUG
#endif
#include <cassert>

using namespace std;
using namespace kim;

ostream &operator<<(ostream &os, const VariantNodesIndex::VariantNode &node) {
  os << "'" << node.variant << "' [" << node.in_degree << "]";
  return os;
}

void variant_infos(const VariantNodesIndex::VariantNode &node, const string &k, uint16_t v) {
  cout << node << " (expecting '" << k << "' [" << v << "])" << endl;
  assert(node.variant == k);
  assert(node.in_degree == v);
}

void index_infos(const VariantNodesIndex &idx, bool empty, size_t size) {
  cout << "Variant node index informations:" << endl;
  cout << "- empty status: " << idx.empty()
       << " (expecting " << empty << ")" << endl;
  assert(idx.empty() == empty);
  cout << "- size: " << idx.size()
       << " (expecting " << size << ")" << endl;
  assert(idx.size() == size);
  cout << "- content: " << endl;
  for (VariantNodesIndex::const_iterator it = idx.cbegin();
       it != idx.cend();
       ++it) {
    cout << "  - " << VariantNodesIndex::VariantNode(it) << endl;
  }
}


int main() {

  cout << "Testing VariantNodesIndex" << endl;

  VariantNodesIndex idx;
  index_infos(idx, true, 0);

  cout << "Adding new variant 'V1'" << endl;
  idx.addVariantNode("V1");
  index_infos(idx, false, 1);

  cout << "Adding new variant 'V2'" << endl;
  idx += "V2";
  index_infos(idx, false, 2);

  cout << "Adding existing variant 'V1'" << endl;
  idx += "V1";
  index_infos(idx, false, 2);

  cout << "Retrieving variant node associated to 'V1'" << endl;
  variant_infos(idx.getVariantNode("V1"), "V1", 0);

  cout << "Retrieving variant node associated to 'V1'" << endl;
  variant_infos(idx.getVariantNode("V2"), "V2", 0);

  cout << "Trying to retrieve inexistant variant node associated to 'V3'" << endl;
  variant_infos(idx.getVariantNode("V3"), "", 0);

  cout << "Incrementing the incoming degree of node associated to 'V1'" << endl;
  idx.addVariantNode("V1", true);
  index_infos(idx, false, 2);
  variant_infos(idx["V1"], "V1", 1);

  cout << "Incrementing the incoming degree of node associated to 'V1'" << endl;
  idx.addVariantNode("V1", true);
  index_infos(idx, false, 2);
  variant_infos(idx["V1"], "V1", 2);

  cout << "That's All, Folk!!!" << endl;
  return 0;

}
