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

#include "kmer_variant_edges_subindex.h"

#include "config.h"
#include "sort_helper.h"

using namespace std;
using namespace bm;

BEGIN_KIM_NAMESPACE

#define ERROR_MSG(msg)                          \
  do {                                          \
    KmerVariantEdgesSubindexException error;    \
    error << msg;                               \
    throw error;                                \
  } while (0)

#define CHECK_FROZEN_STATE(expected_state, mth)                         \
  if (!(expected_state)) {                                              \
    ERROR_MSG("Edges sub-index must be " << (expected_state ? "frozen" : "unfrozen") \
              << " before calling KmerVariantEdgesSubindex::" << #mth   \
              << "() method.");                                         \
  }                                                                     \
  (void) 0

KmerVariantEdgesSubindex::KmerVariantEdgesSubindex(KmerNodesSubindex &kmer_nodes,
                                                   VariantNodesIndex &variant_nodes,
                                                   size_t estimated_nb_edges):
  _first_kmer_edges(),
  _rank_select_helper(),
  _kmer_variant_edges(),
  _kmer_nodes(kmer_nodes),
  _variant_nodes(variant_nodes),
  _max_associations(0),
  _frozen(false)
{
  _first_kmer_edges.resize(0);
  _kmer_variant_edges.reserve(estimated_nb_edges);
}

void KmerVariantEdgesSubindex::compact() {
  CHECK_FROZEN_STATE(!frozen(), compact);
  if (!_kmer_nodes.sorted()) {
    _kmer_nodes.expand(_first_kmer_edges);
    if (_kmer_nodes.size() != size()) {
      ERROR_MSG("The k-mer nodes sub-index and this edges sub-index must have the same size to be compacted.");
    }
    vector<size_t> permutation = _kmer_nodes.sort();
    SortHelper<_KmerVariantEdge>::sort<_KmerVariantEdge>(_kmer_variant_edges, permutation);
    _first_kmer_edges = _kmer_nodes.unique();
  }
  _updateMaxAssociations();
}

void KmerVariantEdgesSubindex::freeze() {
  assert(_kmer_nodes.sorted());
  if (!_first_kmer_edges.is_ro()) {
    _first_kmer_edges.freeze();
  }
  _first_kmer_edges.build_rs_index(&(_rank_select_helper));
  _kmer_nodes.freeze();
  _frozen = true;
}

void KmerVariantEdgesSubindex::unfreeze() {
  if (!frozen()) return;
  if (_first_kmer_edges.is_ro()) {
    bvector<> tmp(_first_kmer_edges, finalization::READWRITE);
    _first_kmer_edges.swap(tmp);
  }
  _kmer_nodes.unfreeze();
  _frozen = false;
}

void KmerVariantEdgesSubindex::reserve(size_t n) {
  CHECK_FROZEN_STATE(!frozen(), reserve);
  size_t old_n = _first_kmer_edges.size();
  _kmer_variant_edges.reserve(n);
  _first_kmer_edges.resize(_kmer_variant_edges.capacity());
  if (old_n < n) {
    _first_kmer_edges.set_range(old_n, n - 1);
  }
}

void KmerVariantEdgesSubindex::shrink_to_fit() {
  CHECK_FROZEN_STATE(!frozen(), shrink_to_fit);
  _kmer_variant_edges.shrink_to_fit();
  _first_kmer_edges.resize(_kmer_variant_edges.size());
}

KmerVariantEdgesSubindex::Edge KmerVariantEdgesSubindex::operator[](size_t i) const {
  size_t p;
  if (frozen()) {
    p = _first_kmer_edges.rank(i, _rank_select_helper);
  } else {
    p = _first_kmer_edges.count_range_no_check(0, i);
  }
  assert(p);
  --p;
  return {
          _kmer_nodes[p],
          _kmer_variant_edges[i].rank,
          VariantNodesIndex::VariantNode(_kmer_variant_edges[i].variant_node_iterator)
  };
}

list<KmerVariantEdgesSubindex::KmerVariantAssociation> KmerVariantEdgesSubindex::getKmerVariantAssociation(const BoundedSizeString &suffix) const {
  CHECK_FROZEN_STATE(frozen(), getKmerVariantAssociation);
  list<KmerVariantAssociation> l;
  size_t rank = _kmer_nodes.getKmerNodeRank(suffix);
  size_t p = _selectKmerFirstEdgePosition(rank);
  if (p != size_t(-1)) {
    do {
      l.emplace_back(_kmer_variant_edges[p]);
    } while ((++p < _kmer_variant_edges.size()) && !_first_kmer_edges.get_bit(p));
  }
  return l;
}

KmerVariantEdgesSubindex &KmerVariantEdgesSubindex::add(KmerNodesSubindex::KmerNode node,
                                                        const string &variant,
                                                        size_t rank) {
  CHECK_FROZEN_STATE(!frozen(), add);
  bool new_node = _kmer_nodes.add(node);
  VariantNodesIndex::iterator variant_node_iterator = _variant_nodes.addVariantNode(variant, true);
  if (rank > uint64_t(uint16_t(-1))) {
    ERROR_MSG("The current k-mer suffix (" << node.suffix
              << ") rank is greater than the maximal size allowed by the current implementation ("
              << (size_t) (uint16_t(-1)) << ")");
  }
  _KmerVariantEdge edge = { variant_node_iterator, uint16_t(rank) };
  _kmer_variant_edges.push_back(edge);
  assert(_kmer_nodes.size() > 0);

  _first_kmer_edges.set_bit(_kmer_variant_edges.size() - 1, new_node);
  return *this;
}

KmerVariantEdgesSubindex &KmerVariantEdgesSubindex::merge(KmerVariantEdgesSubindex &idx) {
  CHECK_FROZEN_STATE(!frozen(), merge);
#ifndef NDEBUG
  size_t orig_n1 = _kmer_nodes.size();
  size_t n2 = idx._kmer_nodes.size();
#endif
  size_t orig_m1 = size();
  assert(orig_m1 >= orig_n1);

  size_t m2 = idx.size();
  assert(m2 >= n2);

  assert(&(idx._variant_nodes) == &_variant_nodes);
  bvector<> pos1, pos2;
  _kmer_nodes.merge(idx._kmer_nodes, pos1, pos2);
  assert(pos1.size() == pos2.size());

  size_t n1 = _kmer_nodes.size();
  assert(pos1.size() == n1);
  assert(n1 >= orig_n1);
  assert(n1 >= n2);
  assert(n1 <= orig_n1 + n2);

  size_t m1 = orig_m1 + m2;
  assert(m1 >= n1);

  _kmer_variant_edges.resize(m1);
  _first_kmer_edges.resize(m1);

  size_t i = orig_m1, j = m2, k = n1;

  while (k--) {
    assert(i + j > 0);
    if (pos2[k]) {
      // Copy the newly inserted edges for the pre-existing node.
      do {
        assert(j);
        --j;
        _kmer_variant_edges[i + j] = idx._kmer_variant_edges[j];
        // Increasing the incoming degree of the variant node.
        _kmer_variant_edges[i + j].variant_node_iterator->second++;
        _first_kmer_edges.set_bit_no_check(i + j, 0);
      } while (!idx._first_kmer_edges.get_bit(j));
    }
    if (j && pos1[k]) {
      // If j == 0 there is no need to move edges.
      // Move the pre-existing associated edges.
      do {
        assert(i);
        --i;
        _kmer_variant_edges[i + j] = _kmer_variant_edges[i];
        assert((i > 0) || _first_kmer_edges.get_bit(i));
        _first_kmer_edges.set_bit_no_check(i + j, 0);
      } while (i && !_first_kmer_edges.get_bit(i));
    }
    _first_kmer_edges.set_bit_no_check(i + j, 1);
  }
  assert((i == 0) || (pos1.count_range(0, i) == 1));
  assert(j == 0);
  _updateMaxAssociations();
  return *this;
}

void KmerVariantEdgesSubindex::clear() {
  CHECK_FROZEN_STATE(!frozen(), clear);
  _first_kmer_edges.resize(0);
  _rank_select_helper.resize(0);
  for (vector<_KmerVariantEdge>::iterator it = _kmer_variant_edges.begin();
       it != _kmer_variant_edges.end();
       ++it) {
    assert(it->variant_node_iterator->second);
    if (--(it->variant_node_iterator->second) == 0) {
      _variant_nodes.erase(it->variant_node_iterator);
    }
  }
  _kmer_variant_edges.clear();
  _kmer_variant_edges.shrink_to_fit();
  _kmer_nodes.clear();
  _kmer_nodes.shrink_to_fit();
  _max_associations = 0;
}

static void multilineString2Stream(ostream &os,
                                   const string &s,
                                   const char delim = '\n',
                                   const string prepend = "#",
                                   const string append = "") {
  size_t n = s.size();
  os << prepend;
  for (size_t i = 0; i < n; ++i) {
    if (s[i] == delim) {
      os << append;
    }
    os << s[i];
    if (s[i] == delim) {
      os << prepend;
    }
  }
  os << append;
}

void KmerVariantEdgesSubindex::toStream(ostream &os, const string &header, const string &footer, bool infos_in_header, bool legend_in_header) const {

  size_t kmer_pos = 0, n = size();

  if (!header.empty()) {
    multilineString2Stream(os, header);
    os << endl;
  }
  if (infos_in_header) {
    os << "#NbKmers: " << _kmer_nodes.size() << endl
       << "#NbAssociations: " << size() << endl
       << "#MaxAssociations: " << getMaxAssociations() << endl;
  }
  if (legend_in_header) {
    os << "#FirstColumn[string]: k-mer suffix" << endl
       << "#SecondColumn[string]: associated variant" << endl
       << "#ThirdColumn[integer]: k-mer rank of the variant" << endl
       << "#FourthColumn[integer]: variant node incoming degree" << endl
       << "#FifthColumn[bit]: k-mer belongs to some reference sequence" << endl;
  }
  for (size_t i = 0; i < n; ++i) {
    kmer_pos += _first_kmer_edges.get_bit(i);
    KmerNodesSubindex::KmerNode kmer_node = _kmer_nodes[kmer_pos - 1];
    KmerVariantAssociation edge = _kmer_variant_edges[i];
    os << kmer_node.suffix << "\t"
       << edge.variant_node.variant << "\t"
       << edge.rank << "\t"
       << kmer_node.in_reference << endl;
  }
  if (!header.empty()) {
    multilineString2Stream(os, footer);
    os << endl;
  }
}

void KmerVariantEdgesSubindex::_updateMaxAssociations() {
  CHECK_FROZEN_STATE(!frozen(), _updateMaxAssociations);
  size_t i = 0, prev = 0, n = _kmer_variant_edges.size();
  _max_associations = (n > 0);
  assert((_first_kmer_edges.size() == 0) || (_first_kmer_edges.count() > 0));
  while (i < n) {
    if (i - prev > _max_associations) {
      _max_associations = i - prev;
    }
    prev = i;
    i = _first_kmer_edges.get_next(prev);
    if (!i) {
      i = n;
    }
  }
  if (prev) {
    assert(i > 0);
    if (i - prev > _max_associations + 1) {
      _max_associations = i - prev;
    }
  }
}

size_t KmerVariantEdgesSubindex::_selectKmerFirstEdgePosition(size_t pos) const {
  CHECK_FROZEN_STATE(frozen(), _selectKmerFirstEdgePosition);
  bvector<>::size_type p;
  if (!_first_kmer_edges.select(pos + 1, p, _rank_select_helper)) {
    p = -1;
  }
  return (p == bvector<>::size_type(-1)) ? size_t(-1) : p;
}

END_KIM_NAMESPACE
