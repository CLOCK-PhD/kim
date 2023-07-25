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

#include <kmer_nodes_subindex.h>
#include <bounded_size_string.h>

#include <string>
#include <vector>
#include <iostream>
#ifdef NDEBUG
#  unef NDEBUG
#endif
#include <cassert>

#include <bm.h>

using namespace std;
using namespace bm;
using namespace kim;

void infos(const KmerNodesSubindex &idx, bool frozen, bool sorted, bool empty, size_t size, size_t capacity) {
  cout << "Sub-index informations:" << endl
       << "- frozen status: " << idx.frozen() << endl;
  assert(idx.frozen() == frozen);
  cout << "- sorted status: " << idx.sorted() << endl;
  assert(idx.sorted() == sorted);
  cout << "- empty status: " << idx.empty() << endl;
  assert(idx.empty() == empty);
  cout << "- size: " << idx.size() << endl;
  assert(idx.size() == size);
  cout << "- capacity: " << idx.capacity() << endl;
  assert(idx.capacity() == capacity);
  cout << "- content: #JSON: " << idx << endl;
  for (size_t i = 0; i < idx.size(); ++i) {
    cout << "  - ";
    idx[i].toYaml(cout);
    cout << endl;
  }
}

#define TEST_EXPECTED_EXCEPTION(code) \
  exception_thrown = false;                                             \
  try code catch (const KmerNodesSubindexException &e) {                \
      cout << "A KmerNodesSubindexException with message '" << e.what() \
           << " has been thrown" << endl;                               \
      exception_thrown = true;                                          \
    }                                                                   \
  assert(exception_thrown)


int main(int argc, char **argv) {

  cout << "Testing KmerNodesSubindex for 3-mers" << endl;

  BoundedSizeString::setMaximalSize(3);
  bool exception_thrown;

  cout << "Creation of an empty sub-index with allocated space for 10 nodes" << endl;
  KmerNodesSubindex idx(10);
  infos(idx, false, true, true, 0, 10);

  cout << "Shrinking to fit the right size" << endl;
  idx.shrink_to_fit();
  infos(idx, false, true, true, 0, 0);

  cout << "Reserving 8 nodes" << endl;
  idx.reserve(8);
  infos(idx, false, true, true, 0, 8);

  cout << "Freezing the sub-index" << endl;
  idx.freeze();
  infos(idx, true, true, true, 0, 8);

  TEST_EXPECTED_EXCEPTION({
      cout << "Trying to reserve space for the frozen sub-index" << endl;
      idx.reserve(10);
    });
  infos(idx, true, true, true, 0, 8);

  TEST_EXPECTED_EXCEPTION({
      cout << "Trying to shrink to fit minimal space for the frozen sub-index" << endl;
      idx.shrink_to_fit();
    });
  infos(idx, true, true, true, 0, 8);

  TEST_EXPECTED_EXCEPTION({
      cout << "Trying to add a kmer node to the frozen sub-index" << endl;
      idx.add(KmerNodesSubindex::KmerNode("AAA", false));
    });
  infos(idx, true, true, true, 0, 8);

  cout << "Unfreezing the sub-index" << endl;
  idx.unfreeze();
  infos(idx, false, true, true, 0, 8);

  vector<BoundedSizeString> kmers = {"AAC", "ACA", "ACT", "AGG", "AGG", "ATA", "TTA" };
  cout << "Adding ordered 3-mers (with some duplications):" << endl;
  for (size_t i = 0; i < kmers.size(); ++i) {
    cout << "- Adding '" << kmers[i] << "'" << endl;
    idx += KmerNodesSubindex::KmerNode(kmers[i], false);
    cout << "- Sub-index front element is '" << idx.front() << "'" << endl;
    assert(idx.front().suffix == kmers[0]);
    cout << "- Sub-index back element is '" << idx.back() << "'" << endl;
    assert(idx.back().suffix == kmers[i]);
  }
  infos(idx, false, true, false, 6, 8);

  cout << "Adding the same 3-mers to the sub-index:" << endl;
  for (size_t i = 0; i < kmers.size(); ++i) {
    cout << "- Adding '" << kmers[i] << "'" << endl;
    bool new_kmer = idx.add(KmerNodesSubindex::KmerNode(kmers[i], false));
    assert((i == 4) xor new_kmer);
  }
  cout << "Shrinking the sub-index" << endl;
  idx.shrink_to_fit();
  infos(idx, false, false, false, 12, 12);
  cout << "Sorting the sub-index" << endl;
  idx.sort();
  infos(idx, false, true, false, 12, 12);
  cout << "Removing duplicated 3-mers from the sub-index" << endl;
  idx.unique();
  infos(idx, false, true, false, 6, 12);

  vector<BoundedSizeString> kmers2 = {"TAC", "ACA", "AAA", "GGG" };
  cout << "Adding unordered 3-mers:" << endl;
  for (size_t i = 0; i < kmers2.size(); ++i) {
    cout << "Adding 3-mer '" << kmers2[i] << "'" << endl;
    idx += KmerNodesSubindex::KmerNode(kmers2[i], false);
  }
  cout << "The sub-index capacity (=" << idx.capacity() << ")"
       << " must be at least 10 (size of the sub-index is " << idx.size() << ")"
       << endl;
  assert(idx.size() == 10);
  assert(idx.capacity() >= 12);
  cout << "Shrinking the sub-index" << endl;
  idx.shrink_to_fit();
  infos(idx, false, false, false, 10, 10);

  TEST_EXPECTED_EXCEPTION({
      cout << "Trying to query the sub-index to find rank of the existing 3-mer '" << kmers2[0] << "'" << endl;
      idx.getKmerNodeRank(kmers2[0]);
    });
  infos(idx, false, false, false, 10, 10);

  cout << "Sorting the sub-index" << endl;
  vector<size_t> permutation = idx.sort();
  vector<size_t> expected_permutation = { 8, 0, 1, 7, 2, 3, 4, 9, 6, 5 };
  cout << "The computed permutation must have 10 elements:" << endl;
  assert(permutation.size() == expected_permutation.size());
  for (size_t i = 0; i < permutation.size(); ++i) {
    cout << "permutation[" << i << "] = " << permutation[i]
         << " (expecting " << expected_permutation[i] << ")" << endl;
    assert(permutation[i] == expected_permutation[i]);
  }
  infos(idx, false, true, false, 10, 10);


  cout << "Removing duplicated 3-mers from the sub-index" << endl;
  bvector<> kept = idx.unique();
  bvector<> expected_kept;
  expected_kept.resize(10);
  expected_kept.set();
  expected_kept[3] = false;
  cout << "The computed kept bit vector must have 10 elements:" << endl;
  assert(kept.size() == expected_kept.size());
  for (size_t i = 0; i < kept.size(); ++i) {
    cout << "kept[" << i << "] = " << kept[i]
         << " (expecting " << expected_kept[i] << ")" << endl;
    assert(kept[i] == expected_kept[i]);
  }
  infos(idx, false, true, false, 9, 10);

  cout << "Querying the sub-index to find rank of the existing 3-mer '" << kmers2[0] << "'" << endl;
  size_t pos = idx.getKmerNodeRank(kmers2[0]);
  cout << "=> found rank " << pos << " (expecting 7)" << endl;
  assert(pos == 7);
  infos(idx, false, true, false, 9, 10);

  cout << "Freezing the sub-index" << endl;
  idx.freeze();
  infos(idx, true, true, false, 9, 10);

  cout << "Creating a new k-mer node sub-index." << endl;
  KmerNodesSubindex idx2;
  infos(idx2, false, true, true, 0, 0);
  vector<BoundedSizeString> kmers3 = {"AGT", "CAG", "TGA", "AGT", "TGA", "CCC", "GCG", "GGG" };
  cout << "Adding unordered 3-mers:" << endl;
  for (size_t i = 0; i < kmers3.size(); ++i) {
    cout << "Adding 3-mer '" << kmers3[i] << "'" << endl;
    idx2 += KmerNodesSubindex::KmerNode(kmers3[i], true);
  }
  cout << "The sub-index capacity (=" << idx2.capacity()
       << " must be at least 8 (size of the sub-index is " << idx2.size()
       << ")" << endl;
  assert(idx2.size() == 8);
  assert(idx2.capacity() >= 8);
  idx2.shrink_to_fit();
  infos(idx2, false, false, false, 8, 8);

  cout << "Sorting the sub-index" << endl;
  permutation = idx2.sort();
  expected_permutation = { 0, 3, 1, 5, 6, 7, 2, 4 };
  cout << "The computed permutation must have 8 elements:" << endl;
  assert(permutation.size() == expected_permutation.size());
  for (size_t i = 0; i < permutation.size(); ++i) {
    cout << "permutation[" << i << "] = " << permutation[i]
         << " (expecting " << expected_permutation[i] << ")" << endl;
    assert(permutation[i] == expected_permutation[i]);
  }
  infos(idx2, false, true, false, 8, 8);

  cout << "Removing duplicated 3-mers from the sub-index" << endl;
  kept = idx2.unique();
  expected_kept.resize(8);
  expected_kept.set();
  expected_kept[1] = expected_kept[7] = false;
  cout << "The computed kept bit vector must have 8 elements:" << endl;
  assert(kept.size() == expected_kept.size());
  for (size_t i = 0; i < kept.size(); ++i) {
    cout << "kept[" << i << "] = " << kept[i]
         << " (expecting " << expected_kept[i] << ")" << endl;
    assert(kept[i] == expected_kept[i]);
  }
  infos(idx2, false, true, false, 6, 8);

  // idx =  {"AAA", "AAC", "ACA", "ACT", "AGG", "ATA", "GGG", "TAC", "TTA"}
  // idx2 = {"AGT", "CAG", "CCC", "GCG", "GGG", "TGA" }
  cout << "The 'GGG' k-mer is common to the two sub-indexes but with different in_reference values:" << endl;
  size_t p1 = idx.getKmerNodeRank(BoundedSizeString("GGG"));
  cout << "In first sub-index, this k-mer is at rank " << p1 << " ( expecting 6)" << endl;
  assert(p1 == 6);
  cout << "Its in_reference value is " << idx.isInReferenceKmer(p1) << " (expecting 0)" << endl;
  size_t p2 = idx2.getKmerNodeRank(BoundedSizeString("GGG"));
  cout << "In second sub-index, this k-mer is at rank " << p2 << " ( expecting 4)" << endl;
  assert(p2 == 4);
  cout << "Its in_reference value is " << idx.isInReferenceKmer(p2) << " (expecting 0)" << endl;

  cout << "Unfreezing the sub-index" << endl;
  idx.unfreeze();
  infos(idx, false, true, false, 9, 10);

  cout << "Make a backup of the first k-mer node sub-index." << endl;
  KmerNodesSubindex idx_orig = idx;
  TEST_EXPECTED_EXCEPTION({
      cout << "Trying to merge the two index knowing that at least one common k-mer having in_reference discordant values." << endl;
      idx += idx2;
    });
  cout << "Restoring initial sub-index from backup" << endl;
  idx = idx_orig;
  infos(idx, false, true, false, 9, 10);


  cout << "Setting the in_reference value to 1 for the common k-mer" << endl;
  idx.setInReferenceKmer(BoundedSizeString("GGG"));
  cout << "Now, k-mer node associated to 'GGG' has in_reference value set to " << idx.isInReferenceKmer(p1) << " (expecting 1)" << endl;
  assert(idx.isInReferenceKmer(p1));

  cout << "Merging the two sub-indexes" << endl;
  bvector<> pos1, pos2;
  idx.merge(idx2, pos1, pos2);
  idx.shrink_to_fit();

  // old_idx = {"AAA", "AAC", "ACA", "ACT", "AGG"       , "ATA",                      "GGG", "TAC",        "TTA"}
  // idx2 =    {                                   "AGT",        "CAG", "CCC", "GCG", "GGG",        "TGA" }
  // idx =     {"AAA", "AAC", "ACA", "ACT", "AGG", "AGT", "ATA", "CAG", "CCC", "GCG", "GGG", "TAC", "TGA", "TTA"}

  infos(idx, false, true, false, 14, 14);
  bvector<> expected_pos1; // { 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1};
  expected_pos1.resize(14);
  expected_pos1.set();
  expected_pos1[5] = expected_pos1[7] = expected_pos1[8] = expected_pos1[9] = expected_pos1[12] = false;
  assert(expected_pos1.size() == 14);
  bvector<> expected_pos2; // { 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0};
  expected_pos2 = expected_pos1;
  expected_pos2.flip();
  expected_pos2[10] = true;
  assert(expected_pos2.size() == 14);
  cout << "The position bit vectors must have size 14." << endl;
  assert(pos1.size() == expected_pos1.size());
  assert(pos2.size() == expected_pos2.size());
  for (size_t i = 0; i < pos1.size(); ++i) {
    cout << "pos1[" << i << "] = " << pos1[i]
         << " (expecting " << expected_pos1[i] << ")\t"
         << "pos2[" << i << "] = " << pos2[i]
         << " (expecting " << expected_pos2[i] << ")" << endl;
    assert(pos1[i] == expected_pos1[i]);
    assert(pos2[i] == expected_pos2[i]);
  }
  infos(idx, false, true, false, 14, 14);

  cout << "Clearing the second sub-index" << endl;
  idx2.clear();
  idx2.shrink_to_fit();
  infos(idx2, false, true, true, 0, 0);

  cout << "That's All, Folk!!!" << endl;
  return 0;

}
