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

#include "kmer_nodes_subindex.h"

#include "config.h"
#include "sort_helper.h"

using namespace bm;

BEGIN_KIM_NAMESPACE

#define ERROR_MSG(msg)                          \
  do {                                          \
    KmerNodesSubindexException error;           \
    error << msg;                               \
    throw error;                                \
  } while (0)

#define CHECK_FROZEN_STATE(expected_state, mth)                         \
  if (!(expected_state)) {                                              \
    ERROR_MSG("k-mer nodes sub-index must be " << (expected_state ? "frozen" : "unfrozen") \
              << " before calling KmerNodesSubindex::" << #mth          \
              << "() method.");                                         \
  }                                                                     \
  (void) 0

KmerNodesSubindex::KmerNodesSubindex(size_t estimated_nb_kmers, bool sorted):
  _suffixes(), _in_reference(), _sorted(sorted), _frozen(false)
{
  _suffixes.reserve(estimated_nb_kmers);
  _in_reference.resize(0);
}

void KmerNodesSubindex::freeze() {
  if (frozen()) return;
  assert(_in_reference.size() == _suffixes.size());
  sort();
  _in_reference.freeze();
  _frozen = true;
}

void KmerNodesSubindex::unfreeze() {
  if (!frozen()) return;
  if (_in_reference.is_ro()) {
    bvector<> tmp(_in_reference, finalization::READWRITE);
    _in_reference.swap(tmp);
  }
  _frozen = false;
}

void KmerNodesSubindex::reserve(size_t n) {
  CHECK_FROZEN_STATE(!frozen(), reserve);
  _suffixes.reserve(n);
}

void KmerNodesSubindex::shrink_to_fit() {
  CHECK_FROZEN_STATE(!frozen(), shrink_to_fit);
  _suffixes.shrink_to_fit();
}

bool KmerNodesSubindex::add(const KmerNode &node) {
  CHECK_FROZEN_STATE(!frozen(), add);
  bool to_add = true;
  size_t p = _suffixes.size();
  if (p) {
    --p;
    int cmp = node.suffix.compare(_suffixes[p]);
    if (cmp == 0) {
      to_add = !_sorted;
      if (node.in_reference != _in_reference[p]) {
        ERROR_MSG("Unable to add the k-mer suffix '" << node.suffix
                  << "' to the current sub-index since it already exists and it was said to be "
                  << (_in_reference[p] ? "present" : "absent")
                  << " from reference sequences and now the opposite is said.");
      }
    } else {
      _sorted &= (cmp > 0);
      ++p;
    }
  }
  if (to_add) {
    _suffixes.push_back(node.suffix);
    _in_reference[p] = node.in_reference;
  }
  return to_add;
}

void KmerNodesSubindex::setInReferenceKmer(size_t rank, bool status) {
  CHECK_FROZEN_STATE(!frozen(), setInReferenceKmer);
  if (rank >= _suffixes.size()) {
    ERROR_MSG("Unable to change the status of the k-mer node at rank " << rank
              << " since there is only " << _suffixes.size() << " nodes.");
  }
  _in_reference[rank] = status;
}

bool KmerNodesSubindex::isInReferenceKmer(size_t rank) const {
  if (rank >= _suffixes.size()) {
    ERROR_MSG("Unable to get the status of the k-mer node at rank " << rank
              << " since there is only " << _suffixes.size() << " nodes.");
  }
  return _in_reference[rank];
}

size_t KmerNodesSubindex::getKmerNodeRank(const BoundedSizeString &suffix) const {
  if (!_sorted) {
    ERROR_MSG("Unable to perform a dichotomic lookup in the k-mer node sub-index since it is not sorted");
  }

  size_t start = 0, end = _suffixes.size(), p;
  int cmp = 1;
  while (start < end) {
    p = (start + end) / 2;
    cmp = _suffixes[p].compare(suffix);
    if (cmp < 0) {
      start = p + 1;
    } else {
      if (cmp > 0) {
        end = p;
      } else {
        start = end = p;
      }
    }
  }
  return cmp ? size_t(-1) : p;
}

std::vector<size_t> KmerNodesSubindex::sort() {
  CHECK_FROZEN_STATE(!frozen(), sort);
  SortHelper<BoundedSizeString> helper(_suffixes);
  helper.sort<bool, bvector<>, bvector<>::reference>(_in_reference);
  helper.sort<BoundedSizeString>(_suffixes);
  _sorted = true;
  return helper.permutation();
}

bvector<> KmerNodesSubindex::unique() {
  CHECK_FROZEN_STATE(!frozen(), unique);
  size_t i = 0, j = 1, n = _suffixes.size();
  bvector<> kept;
  kept.resize(n);
  kept.set();
  while (j < n) {
    if (_suffixes[i].compare(_suffixes[j])) {
      if (++i < j) {
        _suffixes[i] = _suffixes[j];
        _in_reference[i] = _in_reference[j];
      }
    } else {
      // suffixes at rank i an j are the same.
      if (_in_reference[i] != _in_reference[j]) {
        ERROR_MSG("Unable to remove duplicated k-mer suffix '" << _suffixes[j]
                  << "' from the current sub-index since it already exists and it was said to be "
                  << (_in_reference[i] ? "present" : "absent")
                  << " from reference sequences and now the opposite is said.");
      }
      kept[j] = false;
    }
    ++j;
  }
  if (j > ++i) {
    _suffixes.resize(i);
    _in_reference.resize(i);
  }
  return kept;
}

void KmerNodesSubindex::expand(const bvector<> &schema) {
  size_t orig_n = _suffixes.size();
  assert(_in_reference.size() == orig_n);
  if (orig_n == 0) return;
  size_t n = schema.size();
  assert(n >= orig_n);
  assert(schema.count() == orig_n);
  _suffixes.resize(n);
  _in_reference.resize(n);
  --orig_n;
  while (orig_n < --n) {
    _suffixes[n] = _suffixes[orig_n];
    _in_reference[n] = _in_reference[orig_n];
    if (schema[n]) {
      --orig_n;
    }
  }
}

KmerNodesSubindex &KmerNodesSubindex::merge(const KmerNodesSubindex &idx, bvector<> &pos1, bvector<> &pos2) {
  CHECK_FROZEN_STATE(!frozen(), merge);
  if (!sorted()) {
    ERROR_MSG("Current sub-index must be sorted before calling the KmerNodesSubindex::merge() method.");
  }
  if (!idx.sorted()) {
    ERROR_MSG("Right hand side sub-index must be sorted before calling the KmerNodesSubindex::merge() method.");
  }
  // First pass, compute the number of common nodes.
  size_t i = 0, j = 0, k = 0, m = size(), n = idx.size();
  pos1.clear();
  pos2.clear();
  pos1.resize(m+n);
  pos2.resize(m+n);
  while ((i < m) && (j < n)) {
    int cmp = _suffixes[i].compare(idx._suffixes[j]);
    if (cmp < 0) {
      pos1[k] = true;
      ++i;
    } else {
      if (cmp == 0) {
        if (_in_reference[i] != idx._in_reference[j]) {
          ERROR_MSG("Unable to add the k-mer suffix '" << idx._suffixes[j]
                    << "' to the current sub-index since it already exists and it was said to be "
                    << (_in_reference[i] ? "present" : "absent")
                    << " from reference sequences and now the opposite is said.");
        }
        pos1[k] = true;
        ++i;
      }
      pos2[k] = true;
      ++j;
    }
    ++k;
  }

  if (i < m) {
    pos1.set_range_no_check(k, k + m - i - 1);
    k += m - i;
    i = m;
  }
  if (j < n) {
    pos2.set_range_no_check(k, k + n - i - 1);
    k += n - j;
    i = n;
  }
  pos1.resize(k);
  pos2.resize(k);
  _suffixes.resize(k);
  _in_reference.resize(k);
  _suffixes.shrink_to_fit();

  // Second pass, merge idx values into current sub-index starting by
  assert(i == m);
  assert(j == n);
  assert(k <= i + j);
  while (i && j) {
    assert(k > 0);
    if (pos1[--k]) {
      // The existing k-mer node at rank (i - 1) in this sub-index
      // must be moved to rank k
      _suffixes[k] = _suffixes[--i];
      _in_reference[k] = _in_reference[i];
      if (pos2[k]) {
        --j;
      }
    } else {
      // Whatever the situation, j is decreased by onea dn pos2[j] is
      // true by construction.
      --j;
      assert(pos2[k]);
      // The k-mer at rank j in the sub-index to merge must be added
      // at rank k in this sub-index.
      _suffixes[k] = idx._suffixes[j];
      _in_reference[k] = idx._in_reference[j];
        // k is decreased by one. A false value for 'kept[--j]' in the
        // conditional means that the k-mer at rank j in the index to merge is
        // already is in this sub-index and must simply been ignored.
    }
  }

  if (j--) {
    --k;
    assert(j == k);
    // All k-mers until rank j in the index to merge must be moved in front of the current
    move(idx._suffixes.begin(), idx._suffixes.begin() + i, _suffixes.begin());
    _in_reference.copy_range(idx._in_reference, 0, j);
    j -= k;
    // A zero value for j in the conditional means that
    // all k-mers until rank i are already sorted and at their final
    // ranks (there is nothing to do).
  }
  assert(k == i);
  return *this;
}

void KmerNodesSubindex::clear() {
  CHECK_FROZEN_STATE(!frozen(), clear);
  _suffixes.clear();
  _in_reference.resize(0);
  _sorted = true;
}

END_KIM_NAMESPACE
