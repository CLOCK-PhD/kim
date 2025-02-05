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

#include <kmer_variant_edges_subindex.h>

#include <iostream>
#ifdef NDEBUG
#  undef NDEBUG
#endif
#include <cassert>

#include <bm.h>

using namespace std;
using namespace bm;
using namespace kim;

ostream &operator<<(ostream &os, const KmerVariantEdgesSubindex::Edge &e) {
  os << e.kmer_node.suffix;
  if (e.kmer_node.in_reference) {
    os << "[ref]";
  }
  os << " to " << e.variant_node.variant
     << " (in_degree: " << e.variant_node.in_degree
     << ", rank: " << e.rank << ")";
  return os;
}

void infos(const KmerVariantEdgesSubindex &idx, bool frozen, bool empty, size_t size, size_t capacity, size_t max_associations, size_t nb_kmer_nodes, size_t nb_variant_nodes) {
  cout << "Sub-index informations:" << endl
       << "- frozen status: " << idx.frozen()
       << " (expecting " << frozen << ")" << endl;
  assert(idx.frozen() == frozen);
  cout << "- empty status: " << idx.empty()
       << " (expecting " << empty << ")" << endl;
  assert(idx.empty() == empty);
  cout << "- size: " << idx.size()
       << " (expecting " << size << ")" << endl;
  assert(idx.size() == size);
  cout << "- capacity: " << idx.capacity()
       << " (expecting " << capacity << ")" << endl;
  assert(idx.capacity() == capacity);
  cout << "- maximal number of associations: " << idx.getMaxAssociations()
       << " (expecting " << max_associations << ")" << endl;
  assert(idx.getMaxAssociations() == max_associations);
  cout << "- nb k-mer nodes: " << idx.kmerNodesSubindex().size()
       << " (expecting " << nb_kmer_nodes << ")" << endl;
  assert(idx.kmerNodesSubindex().size() == nb_kmer_nodes);
  cout << "- nb variant nodes: " << idx.variantNodesIndex().size()
       << " (expecting " << nb_variant_nodes << ")" << endl;
  assert(idx.variantNodesIndex().size() == nb_variant_nodes);
  cout << "- content of edges associations: ";
  for (size_t i = 0; i < idx.size(); ++i) {
    cout << "  - " << idx[i] << endl;
  }
  cout << "- content of sub-index: " << endl;
  idx.toStream(cout);
  cout << endl;
}

#define TEST_EXPECTED_EXCEPTION(code)                                   \
  exception_thrown = false;                                             \
  try code catch (const KmerVariantEdgesSubindexException &e) {         \
      cout << "A KmerVariantEdgesSubindexException with message '" << e.what() \
           << "' has been thrown" << endl;                              \
      exception_thrown = true;                                          \
    }                                                                   \
    catch (const KmerNodesSubindexException &e) {                       \
      cout << "A KmerNodesSubindexException with message '" << e.what() \
           << "' has been thrown" << endl;                              \
      exception_thrown = true;                                          \
    }                                                                   \
  cout << endl;                                                         \
  assert(exception_thrown)


int main() {

  cout << "Testing KmerVariantEdgesSubindex" << endl;

  bool exception_thrown;
  BoundedSizeString::setMaximalSize(3);
  KmerNodesSubindex kmer_nodes;
  VariantNodesIndex variant_nodes;

  cout << "Creation of an empty sub-index expecting edges to be added into"
    " lexicographic order of the k-mer source node (no k-mer nodes and no"
    " variant nodes at beginning and space for 10 edges)" << endl;
  KmerVariantEdgesSubindex idx(kmer_nodes, variant_nodes, 10);
  infos(idx, false, true, 0, 10, 0, 0, 0);

  cout << "Shrinking to fit the right size" << endl;
  idx.shrink_to_fit();
  infos(idx, false, true, 0, 0, 0, 0, 0);

  cout << "Reserving 8 edges" << endl;
  idx.reserve(8);
  infos(idx, false, true, 0, 8, 0, 0, 0);

  cout << "Freezing the sub-index" << endl;
  idx.freeze();
  infos(idx, true, true, 0, 8, 0, 0, 0);

  TEST_EXPECTED_EXCEPTION({
      cout << "Trying to reserve space for the frozen sub-index" << endl;
      idx.reserve(10);
    });
  infos(idx, true, true, 0, 8, 0, 0, 0);

  TEST_EXPECTED_EXCEPTION({
      cout << "Trying to shrink to fit minimal space for the frozen sub-index" << endl;
      idx.shrink_to_fit();
    });
  infos(idx, true, true, 0, 8, 0, 0, 0);

  TEST_EXPECTED_EXCEPTION({
      cout << "Trying to add an edge to the frozen sub-index" << endl;
      idx.add(KmerNodesSubindex::KmerNode("AAA", false), "Var1", 0);
    });
  infos(idx, true, true, 0, 8, 0, 0, 0);

  cout << "Unfreezing the sub-index" << endl;
  idx.unfreeze();
  infos(idx, false, true, 0, 8, 0, 0, 0);

  cout << "Adding an edge from new k-mer node 'AAA' to new variant 'V1' with rank 1" << endl;
  idx.add(KmerNodesSubindex::KmerNode("AAA", false), "V1", 1);
  /*
   * Expected sub-index idx is:
   * k-mers: AAA
   * variants: V1[1]
   * edges:
   * - AAA -> V1 [1]
   */
  infos(idx, false, false, 1, 8, 0, 1, 1);

  cout << "Adding an edge from existing k-mer node 'AAA' to new variant 'V2' with rank 0" << endl;
  idx.add(KmerNodesSubindex::KmerNode("AAA", false), "V2", 0);
  /*
   * Expected sub-index idx is:
   * k-mers: AAA
   * variants: V1[1] V2[1]
   * edges:
   * - AAA -> V1 [1]
   * - AAA -> V2 [0]
   */
  infos(idx, false, false, 2, 8, 0, 1, 2);

  cout << "Adding an edge from new k-mer node 'CCC' to existing variant 'V1' with rank 2" << endl;
  idx.add(KmerNodesSubindex::KmerNode("CCC", false), "V1", 2);
  /*
   * Expected sub-index idx is:
   * k-mers: AAA CCC
   * variants: V1[2] V2[1]
   * edges:
   * - AAA -> V1 [1]
   * - AAA -> V2 [0]
   * - CCC -> V1 [2]
   */
  infos(idx, false, false, 3, 8, 0, 2, 2);

  cout << "Adding an edge from new k-mer node 'GGG' to new variant 'V3' with rank 1" << endl;
  idx.add(KmerNodesSubindex::KmerNode("GGG", false), "V3", 1);
  /*
   * Expected sub-index idx is:
   * k-mers: AAA CCC GGG
   * variants: V1[2] V2[1] V3[1]
   * edges:
   * - AAA -> V1 [1]
   * - AAA -> V2 [0]
   * - CCC -> V1 [2]
   * - GGG -> V3 [1]
   */
  infos(idx, false, false, 4, 8, 0, 3, 3);

  cout << "Adding an edge from new k-mer node 'TTT' to existing variant 'V1' with rank 4" << endl;
  idx.add(KmerNodesSubindex::KmerNode("TTT", false), "V1", 4);
  /*
   * Expected sub-index idx is:
   * k-mers: AAA CCC GGG TTT
   * variants: V1[3] V2[1] V3[1]
   * edges:
   * - AAA -> V1 [1]
   * - AAA -> V2 [0]
   * - CCC -> V1 [2]
   * - GGG -> V3 [1]
   * - TTT -> V1 [4]
   */
  infos(idx, false, false, 5, 8, 0, 4, 3);

  cout << "Compacting the sub-index" << endl;
  idx.compact();
  /*
   * Expected sub-index idx is:
   * k-mers: AAA CCC GGG TTT
   * variants: V1[3] V2[1] V3[1]
   * edges:
   * - AAA -> V1 [1]
   * - AAA -> V2 [0]
   * - CCC -> V1 [2]
   * - GGG -> V3 [1]
   * - TTT -> V1 [4]
   */
  infos(idx, false, false, 5, 8, 2, 4, 3);


  cout << "Creation of another sub-index using the same index for variant nodes"
    " as the first edges sub-index but where edges are not expected to be"
    " added in lexicographic order of the k-mer suffix node." << endl;
  KmerNodesSubindex kmer_nodes2;
  KmerVariantEdgesSubindex idx2(kmer_nodes2, variant_nodes, 0);
  /*
   * Expected sub-index idx2 is:
   * k-mers:
   * variants: V1[3] V2[1] V3[1]
   * edges:
   */
  infos(idx2, false, true, 0, 0, 0, 0, 3);

  cout << "Adding an edge from existing k-mer node 'CAA' to existing variant 'V3' with rank 2" << endl;
  idx2.add(KmerNodesSubindex::KmerNode("CAA", false), "V3", 1);
  /*
   * Expected sub-index idx2 is:
   * k-mers: CAA
   * variants: V1[3+0] V2[1+0] V3[1+1]
   * edges:
   * - CAA -> V3 [1]
   */
  infos(idx2, false, false, 1, 1, 0, 1, 3);

  cout << "Adding an edge from existing k-mer node 'ACA' to existing variant 'V4' with rank 5" << endl;
  idx2.add(KmerNodesSubindex::KmerNode("ACA", false), "V4", 5);
  /*
   * Expected sub-index idx2 is:
   * k-mers: CAA ACA
   * variants: V1[3+0] V2[1+0] V3[1+1] V4[0+1]
   * edges:
   * - CAA -> V3 [1]
   * - ACA -> V4 [5]
   */
  infos(idx2, false, false, 2, 2, 0, 2, 4);

  cout << "Adding an edge from existing k-mer node 'ACA' to new variant 'V5' with rank 6" << endl;
  idx2.add(KmerNodesSubindex::KmerNode("ACA", false), "V5", 6);
  /*
   * Expected sub-index idx2 is:
   * k-mers: CAA ACA ACA
   * variants: V1[3+0] V2[1+0] V3[1+1] V4[0+1] V5[0+1]
   * edges:
   * - CAA -> V3 [1]
   * - ACA -> V4 [5]
   * - ACA -> V5 [6]
   */
  infos(idx2, false, false, 3, 4, 0, 3, 5);

  cout << "Adding an edge from existing k-mer node 'CCC' to existing variant 'V4' with rank 8" << endl;
  idx2.add(KmerNodesSubindex::KmerNode("CCC", false), "V4", 8);
  /*
   * Expected sub-index idx2 is:
   * k-mers: CAA ACA ACA CCC
   * variants: V1[3+0] V2[1+0] V3[1+1] V4[0+2] V5[0+1]
   * edges:
   * - CAA -> V3 [1]
   * - ACA -> V4 [5]
   * - ACA -> V5 [6]
   * - CCC -> V4 [8]
   */
  infos(idx2, false, false, 4, 4, 0, 4, 5);

  cout << "Adding an edge from existing k-mer node 'CAA' to existing variant 'V2' with rank 7" << endl;
  idx2.add(KmerNodesSubindex::KmerNode("CAA", false), "V2", 7);
  /*
   * Expected sub-index idx2 is:
   * variants: V1[3+0] V2[1+1] V3[1+1] V4[0+2] V5[0+1]
   * variants: V1 V2 V3 V4 V5
   * edges:
   * - CAA -> V3 [1]
   * - ACA -> V4 [5]
   * - ACA -> V5 [6]
   * - CCC -> V4 [8]
   * - CAA -> V2 [7]
   */
  infos(idx2, false, false, 5, 8, 0, 5, 5);

  cout << "Adding an edge from existing k-mer node 'ACA' to new variant 'V2' with rank 1" << endl;
  idx2.add(KmerNodesSubindex::KmerNode("ACA", false), "V2", 1);
  /*
   * Expected sub-index idx2 is:
   * k-mers: CAA ACA ACA CCC CAA ACA
   * variants: V1[3+0] V2[1+2] V3[1+1] V4[0+2] V5[0+1]
   * edges:
   * - CAA -> V3 [1]
   * - ACA -> V4 [5]
   * - ACA -> V5 [6]
   * - CCC -> V4 [8]
   * - CAA -> V2 [7]
   * - ACA -> V2 [1]
   */
  infos(idx2, false, false, 6, 8, 0, 6, 5);

  cout << "Adding an edge from existing k-mer node 'TTT' to existing variant 'V4' with rank 9" << endl;
  idx2.add(KmerNodesSubindex::KmerNode("TTT", false), "V4", 9);
  /*
   * Expected sub-index idx2 is:
   * k-mers: CAA ACA ACA CCC CAA ACA TTT
   * variants: V1[3+0] V2[1+2] V3[1+1] V4[0+3] V5[0+1]
   * edges:
   * - CAA -> V3 [1]
   * - ACA -> V4 [5]
   * - ACA -> V5 [6]
   * - CCC -> V4 [8]
   * - CAA -> V2 [7]
   * - ACA -> V2 [1]
   * - TTT -> V4 [9]
   */
  infos(idx2, false, false, 7, 8, 0, 7, 5);

  TEST_EXPECTED_EXCEPTION({
      cout << "Trying to merge not compacted sub-indexes" << endl;
      idx += idx2;
    });

  cout << "Compacting the second sub-index" << endl;
  idx2.compact();
  /*
   * Expected sub-index idx2 is:
   * k-mers: ACA CAA CCC TTT
   * variants: V1[3+0] V2[1+2] V3[1+1] V4[0+3] V5[0+1]
   * edges:
   * - ACA -> V4 [5]
   * - ACA -> V5 [6]
   * - ACA -> V2 [1]
   * - CAA -> V3 [1]
   * - CAA -> V2 [7]
   * - CCC -> V4 [8]
   * - TTT -> V4 [9]
   */
  infos(idx2, false, false, 7, 8, 3, 4, 5);

  cout << "Merging the two compacted sub-indexes" << endl;
  idx += idx2;
  /*
   * Expected sub-index idx is:
   * k-mers: AAA ACA CAA CCC GGG TTT
   * variants: V1[3+0] V2[3+2] V3[2+1] V4[3+3] V5[1+1]
   * edges:
   * - AAA -> V1 [1]
   * - AAA -> V2 [0]
   * - ACA -> V4 [5]
   * - ACA -> V5 [6]
   * - ACA -> V2 [1]
   * - CAA -> V3 [1]
   * - CAA -> V2 [7]
   * - CCC -> V1 [2]
   * - CCC -> V4 [8]
   * - GGG -> V3 [1]
   * - TTT -> V1 [4]
   * - TTT -> V4 [9]
   */
  infos(idx, false, false, 12, 12, 3, 6, 5);

  cout << "Clearing second sub-index" << endl;
  idx2.clear();
  infos(idx2, false, true, 0, 0, 0, 0, 5);

  cout << "Adding an edge from new k-mer node 'CAC' to new variant 'V6' with rank 0" << endl;
  idx2.add(KmerNodesSubindex::KmerNode("CAC", false), "V6", 0);
  /*
   * Expected sub-index idx2 is:
   * k-mers: CAC
   * variants: V1[3] V2[3] V3[2] V4[3] V5[1] V6[1]
   * edges:
   * - CAC -> V6 [0]
   */
  infos(idx2, false, false, 1, 1, 0, 1, 6);

  cout << "Clearing again this second sub-index" << endl;
  idx2.clear();
  infos(idx2, false, true, 0, 0, 0, 0, 5);

  cout << "Final index (after freezing)." << endl;
  idx.freeze();
  /*
   * Expected sub-index idx is:
   * k-mers: AAA ACA CAA CCC GGG TTT
   * variants: V1[3] V2[3] V3[2] V4[3] V5[1]
   * edges:
   * - AAA -> V1 [1]
   * - AAA -> V2 [0]
   * - ACA -> V4 [5]
   * - ACA -> V5 [6]
   * - ACA -> V2 [1]
   * - CAA -> V3 [1]
   * - CAA -> V2 [7]
   * - CCC -> V1 [2]
   * - CCC -> V4 [8]
   * - GGG -> V3 [1]
   * - TTT -> V1 [4]
   * - TTT -> V4 [9]
   */
  infos(idx, true, false, 12, 12, 3, 6, 5);

  cout << "That's All, Folk!!!" << endl;
  return 0;

}
