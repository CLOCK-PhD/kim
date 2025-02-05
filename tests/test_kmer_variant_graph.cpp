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

#include <kmer_variant_graph.h>

#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <unistd.h>
#ifdef NDEBUG
#  undef NDEBUG
#endif
#include <cassert>

#include <bm.h>

using namespace std;
using namespace bm;
using namespace kim;

void check_graph_properties(const KmerVariantGraph &g, size_t nb_kmers, size_t nb_variants, size_t nb_edges, bool frozen) {
  const Settings &s = g.settings();
  cout << "Graph properties:" << endl
       << "- Frozen state: " << g.frozen()
       << " (expecting " << frozen << ")" << endl;
  assert(g.frozen() == frozen);
  cout << "- nb of k-mer nodes: " << g.getNbKmers()
       << " (expecting " << (g.frozen() ? "" : " at least ") << nb_kmers << ")" << endl;
  assert((g.getNbKmers() == nb_kmers) || (!g.frozen() && (g.getNbKmers() >= nb_kmers)));
  cout << "- nb of variant nodes: " << g.getNbVariants()
       << " (expecting " << nb_variants << ")" << endl;
  assert(g.getNbVariants() == nb_variants);
  cout << "- nb of edges from k-mer nodes to variant nodes: " << g.getNbKmerVariantEdges()
       << " (expecting " << nb_edges << ")" << endl;
  assert(g.getNbKmerVariantEdges() == nb_edges);
  const KmerVariantGraph::KmerVariantEdgesIndex &edges = g.getKmerVariantEdgesIndex();
  assert(!edges.empty());
  cout << "- nb of sub-indexes: " << edges.size()
       << " (expecting " << (4u << s.getKmerPrefixLength()) << ")" << endl;
  assert(edges.size() == (4u << s.getKmerPrefixLength()));
  cout << endl;
}

// Don't use bvector instead of vector<bool> since intializer constructor doesn't work.
void check_subindex(const KmerVariantGraph &g, const string &p, const vector<KmerVariantGraph::Edge> e, const vector<bool> &in_ref) {
  const Settings &s = g.settings();
  assert(s.getKmerSuffixLength() == BoundedSizeString::getMaximalSize());
  assert(p.size() == s.getKmerPrefixLength());
  size_t v = encode(p, s.getKmerPrefixLength());
  cout << "Checking sub-index #" << v << " associated to prefix '" << p << "'" << endl;
  const KmerVariantEdgesSubindex &edges = g.getKmerVariantEdgesIndex()[v];
  cout << "- Nb of edges: " << edges.size()
       << " (expecting " << e.size() << ")" << endl;
  assert(edges.size() == e.size());
  for (size_t i = 0; i < e.size(); ++i) {
    KmerVariantGraph::Edge a = e[i];
    KmerVariantEdgesSubindex::Edge b = edges[i];
    cout << "  - "
         << a.kmer << "[" << in_ref[i] << "]"
         << " to[" << a.rank  << "] "
         << a.variant
         << " (expecting "
         << b.kmer_node.suffix << "[" << b.kmer_node.in_reference << "]"
         << " to[" << b.rank << "] "
         << b.variant_node.variant << ")" << endl;
    assert(a.kmer.size() == BoundedSizeString::getMaximalSize());
    assert(BoundedSizeString(a.kmer) == b.kmer_node.suffix);
    assert(in_ref[i] == b.kmer_node.in_reference);
    assert(a.rank == b.rank);
    assert(a.variant == b.variant_node.variant);
  }
  cout << endl;
}

void display(const KmerVariantGraph &g) {
  const Settings &s = g.settings();
  cout << "Indexed Variants:";
  const KmerVariantGraph::KmerVariantEdgesIndex &all_edges = g.getKmerVariantEdgesIndex();
  assert(all_edges.size() == g.getNbSubindexes());
  const VariantNodesIndex &variants = g.getVariantNodesIndex();
  for (VariantNodesIndex::const_iterator it = variants.cbegin();
       it != variants.cend();
       ++it) {
    VariantNodesIndex::VariantNode v(it);
    cout << " " << v.variant << "[in_degree:" << v.in_degree << "]";
  }
  cout << endl;
  for (size_t i = 0; i < all_edges.size(); ++i) {
    cout << "Sub-index for k-mer starting by prefix '" << decode(i, s.getKmerPrefixLength()) << "':" << endl;
    const KmerVariantEdgesSubindex &edges = all_edges[i];
    const KmerNodesSubindex &kmer_nodes = edges.kmerNodesSubindex();
    cout << "- (suffix) k-mer nodes: " << kmer_nodes << endl;
    cout << "- Edges:" << endl;
    edges.toStream(cout, "BEGIN", "End", true, false);
  }
  cout << "=====================" << endl << endl;
}

void check_final_graph(const KmerVariantGraph &g, bool frozen) {
  check_graph_properties(g, 14, 4, 15, frozen);
  display(g);
  check_subindex(g, "AA", { }, { });
  check_subindex(g, "AC", { }, { });
  check_subindex(g, "AG", { { "CAG", 0, "V4" }, { "CAT", 4, "V2" } }, {0, 1});
  check_subindex(g, "AT", { { "CAG", 1, "V2" } }, {1});
  check_subindex(g, "CA", { { "CAG", 4, "V3" }, { "GCA", 3, "V2" } }, {0, 1});
  check_subindex(g, "CC", { { "ACA", 3, "V3" } }, {0});
  check_subindex(g, "CG", { { "TTC", 0, "V1" } }, {1});
  check_subindex(g, "CT", { }, { });
  check_subindex(g, "GA", { { "TCA", 0, "V2" } }, {0});
  check_subindex(g, "GC", { { "AGT", 1, "V4" } }, {0});
  check_subindex(g, "GG", { }, { });
  check_subindex(g, "GT", { { "TCC", 1, "V1" }, { "TCC", 0, "V3" } }, {0, 0});
  check_subindex(g, "TA", { }, { });
  check_subindex(g, "TC", { { "AGC", 2, "V2" }, { "CAC", 2, "V3" } }, {1, 0});
  check_subindex(g, "TG", { }, { });
  check_subindex(g, "TT", { { "CCA", 1, "V3" }, { "CCG", 2, "V1" } }, {0, 0});
}

#define TEST_EXPECTED_EXCEPTION(code)                                   \
  exception_thrown = false;                                             \
  try code catch (const KmerVariantGraphException &e) {                 \
      cout << "A KmerVariantGraphException with message '" << e.what()  \
           << "' has been thrown" << endl;                              \
      exception_thrown = true;                                          \
    }                                                                   \
    catch (const KmerVariantGraphParseError &e) {                       \
      cout << "A KmerVariantGraphParseError with message '" << e.what() \
           << "' has been thrown" << endl;                              \
      exception_thrown = true;                                          \
    }                                                                   \
  cout << endl;                                                         \
  assert(exception_thrown)

struct VariantInfo {
  size_t position;
  char variation;
};

int main() {

  bool exception_thrown;
  cout << "Testing KmerVariantGraph" << endl;

  // Graph will index 5-mers with prefixes of length 2 (thus suffixes of length 3).
  // Warning are enabled and settings are not frozen.
  Settings s(5, 2);
  assert(s.valid());
  assert(!s.frozen());

  // Create an empty graph.
  KmerVariantGraph graph(s);

  // Settings should be frozen after graph initialisation.
  assert(s.frozen());

  check_graph_properties(graph, 0, 0, 0, false);


  string ref = "CGATCCGCATCATCAGCCTAGCGTTCTACAGCTCAGCATT";
  // variants:  --T--A--------------------C-----------G-
  //            0         1         2         3         
  //            0123456789012345678901234567890123456789
  vector<VariantInfo> variants = { { 2, 'T' }, { 5, 'A'}, { 26, 'C' }, { 38, 'G' } };

  // reference: CGATCCGCATCATCAGCCTAGCGTTCTACAGCTCAGCATT
  // V1:        --T-------------------------------------
  //            CGTTC                ~~~~~
  //             GTTCC (also associated to V3)
  //              TTCCG
  //
  // reference: CGATCCGCATCATCAGCCTAGCGTTCTACAGCTCAGCATT
  // V2:        -----A----------------------------------
  //             GATCA
  //              ATCAG    ~~~~~
  //               TCAGC    ~~~~~               ~~~~~
  //                CAGCA                        ~~~~~
  //                 AGCAT                        ~~~~~
  //
  // reference: CGATCCGCATCATCAGCCTAGCGTTCTACAGCTCAGCATT
  // V3:        --------------------------C-------------
  //                                  GTTCC (also associated to V1)
  //                                   TTCCA
  //                                    TCCAC
  //                                     CCACA
  //                                      CACAG
  //
  // reference: CGATCCGCATCATCAGCCTAGCGTTCTACAGCTCAGCATT
  // V4:        --------------------------------------G-
  //                                              AGCAG
  //                                               GCAGT
  //

  size_t n = ref.size();
  cout << "                   ";
  for (size_t i = 0; i < n; ++i) {
    cout << ((i % 10) ? "" : " ") << (i / 10);
  }
  cout << endl;
  cout << "                   ";
  for (size_t i = 0; i < n; ++i) {
    cout << ((i % 10) ? "" : " ") << (i % 10);
  }
  cout << endl;
  cout << "Reference sequence:";
  for (size_t i = 0; i < n; ++i) {
    cout << ((i % 10) ? "" : " ") << ref[i];
  }
  cout << endl;
  cout << "Variations:        ";
  size_t v = 0;
  for (size_t i = 0; i < n; ++i) {
    char c = '-';
    if (variants[v].position == i) {
      c = variants[v].variation;
      ++v;
    }
    cout << ((i % 10) ? "" : " ") << c;
  }
  cout << endl;
  size_t cpt = 0;
  for (size_t i = 0; i < variants.size(); ++i) {
    ostringstream name;
    name << "V" << (i + 1);
    size_t p = variants[i].position;
    assert(p < n);
    cout << "Insert k-mers associated to variant '" << name.str()
         << "' at position " << p
         << " (" << ref[p] << " <-> " << variants[i].variation  << ")"
         << endl;
    size_t pos = s.k() - 1;
    if (p < s.k()) {
      pos = p;
      p = 0;
    } else {
      p -= s.k() - 1;
    }
    size_t lg = pos + s.k();
    string superkmer = ref.substr(p, lg);
    assert(superkmer.size() >= s.k());
    superkmer[pos] = variants[i].variation;
    for (size_t j = 0; j < superkmer.size() - s.k() + 1; ++j) {
      string kmer = superkmer.substr(j, s.k());
      cout << "Adding an edge from new k-mer '" << kmer
           << "' to variant '" << name.str()
           << "' with rank " << j
           << endl;
      graph += { kmer, j, name.str() };
      ++cpt;
      check_graph_properties(graph, cpt, (i + 1), cpt, false);
    }
  }

  cout << "Freeze the graph (this will sort and make unique each k-mer nodes in sub-indexes)." << endl;
  graph.freeze();
  check_graph_properties(graph, cpt - 1, variants.size(), cpt, true);

  cout << "Trying to add a new edge to the frozen graph" << endl;
  TEST_EXPECTED_EXCEPTION({
      graph += KmerVariantGraph::Edge({ "AAAAA", 0, "V7" });
    });
  check_graph_properties(graph, cpt - 1, variants.size(), cpt, true);

  cout << "Trying to set in reference k-mers." << endl;
  TEST_EXPECTED_EXCEPTION({
      for (size_t i = 0; i < ref.size() - s.k() + 1; ++i) {
        string kmer = ref.substr(i, s.k());
        graph.setInReferenceKmer(kmer);
      }
    });
  check_graph_properties(graph, cpt - 1, variants.size(), cpt, true);

  cout << "Unfreeze the graph" << endl;
  graph.unfreeze();
  check_graph_properties(graph, cpt - 1, variants.size(), cpt, false);

  cout << "Setting in reference k-mers." << endl
       << " (only k-mer at position 11, 12, 21, 32, 33 and 34 in the reference"
       << " corresponds to k-mer nodes in the graph)" << endl;
  for (size_t i = 0; i < ref.size() - s.k() + 1; ++i) {
    string kmer = ref.substr(i, s.k());
    bool found = graph.setInReferenceKmer(kmer);
    cout << "- k-mer '" << kmer << "' "
         << " (at position " << i << " in the reference) "
         << (found ? " found and set as reference k-mer " : "not found")
         << endl;
    assert(!found xor ((i == 11) || (i == 12) || (i == 21) || (i == 32) || (i == 33) || (i == 34)));
  }
  cout << endl;

  // Expected final graph:
  // Variants: V1[3] V2[5] V3[5] V4[2]
  // NbNodes: 14
  // NbEdges: 15
  // Non-empty sub-indexes:
  // - prefix 'AG':
  //   - nodes(2): CAG(0) CAT(1)
  //   - edges(1):
  //     - CAG V4 (0)
  //     - CAT V2 (4)
  // - prefix 'AT':
  //   - nodes(1): CAG(1)
  //   - edges(1):
  //     - CAG V2 (1)
  // - prefix 'CA':
  //   - nodes(2): CAG(0) GCA(1)
  //   - edges(2):
  //     - CAG V3 (4)
  //     - GCA V2 (3)
  // - prefix 'CC':
  //   - nodes(1): ACA(0)
  //   - edges(1):
  //     - ACA V3 (3)
  // - prefix 'CG':
  //   - nodes(1): TTC(1)
  //   - edges(1):
  //     - TTC V1 (0)
  // - prefix 'GA':
  //   - nodes(1): TCA(0)
  //   - edges(1):
  //     - TCA V2 (0)
  // - prefix 'GC':
  //   - nodes(1): AGT(0)
  //   - edges(1):
  //     - AGT V4 (1)
  // - prefix 'GT':
  //   - nodes(1): TCC(0)
  //   - edges(2):
  //     - TCC V1 (1)
  //     - TCC V3 (0)
  // - prefix 'TC':
  //   - nodes(2): AGC(1) CAC(0)
  //   - edges(2):
  //     - AGC V2 (2)
  //     - CAC V3 (2)
  // - prefix 'TT':
  //   - nodes(2): CCA(0) CCG(0)
  //   - edges(2):
  //     - CCA V3 (1)
  //     - CCG V1 (2)
  check_final_graph(graph, false);

  const string test_metadata = ("\n\n"
                                "Multiple paragraphs in informations might be reformatted in\n"
                                "order to always start paragraphs by \"- \".\n\n"
                                "Multiple lines inside a parargraph might be reformatted in   \n"
                                "order to always start by two spaces.\n\n\n\n"
                                "- A multiline item must preserved\n"
                                "  -- even if it starts by some dash.  \n  \n  \n  \n"
                                "- Empty lines will be removed."
                                "\n\n\n\n\n");
  const string expected_metadata = ("- Multiple paragraphs in informations might be reformatted in\n"
                                    "  order to always start paragraphs by \"- \".\n"
                                    "- Multiple lines inside a parargraph might be reformatted in\n"
                                    "  order to always start by two spaces.\n"
                                    "- A multiline item must preserved\n"
                                    "  -- even if it starts by some dash.\n"
                                    "- Empty lines will be removed.");
  cout << "before setting extra metadata, graph.extraMetadata() should be empty" << endl;
  assert(graph.extraMetadata().empty());

  cout << "Setting the following extra metadata:" << endl
       << "======" << endl
       << test_metadata << endl
       << "======" << endl;
  graph.extraMetadata(test_metadata);
  cout << "Now, the extra metadata are:" << endl
       << "======" << endl
       << graph.extraMetadata() << endl
       << "======" << endl;
  cout << "Expecting the following extra metadata:" << endl
       << "======" << endl
       << expected_metadata << endl
       << "======" << endl;
  assert(graph.extraMetadata() == expected_metadata);

  cout << "After setting empty extra metadata." << endl;
  graph.extraMetadata("");
  cout << "the extra metadata are:" << endl
       << "======" << endl
       << graph.extraMetadata() << endl
       << "======" << endl;
  assert(graph.extraMetadata().empty());

  cout << "Setting the following extra metadata:" << endl
       << "======" << endl
       << expected_metadata << endl
       << "======" << endl;
  graph.extraMetadata(expected_metadata);
  cout << "Now, the extra metadata are:" << endl
       << "======" << endl
       << graph.extraMetadata() << endl
       << "======" << endl;
  cout << "Expecting the following extra metadata:" << endl
       << "======" << endl
       << expected_metadata << endl
       << "======" << endl;
  assert(graph.extraMetadata() == expected_metadata);

  char top_path[] = "test_kmer_variant_graph.XXXXXX";
  if (mkdtemp(top_path)) {
    string path = top_path;
    path += "/kim_index";
    bool first_pass = false;
    do {
      first_pass = !first_pass;
      cout << "Dumping the graph to '" << path << "' index directory." << endl;
      graph.dump(path);
      // Dumping the graph freezes it.
      assert(graph.frozen());
      cout << endl;

      cout << "Unfreezing the graph." << endl;
      graph.unfreeze();
      assert(!graph.frozen());
      cout << endl;

      cout << "Clearing the graph." << endl;
      graph.clear();
      check_graph_properties(graph, 0, 0, 0, false);
      cout << endl;

      cout << "Ensure that extra metadata are empty." << endl;
      assert(graph.extraMetadata().empty());

      cout << "Loading the graph from the dumped index" << endl;
      graph.load(path);
      check_final_graph(graph, true);
      cout << endl;

      cout << "Checking if additional informations are correctly restored:" << endl
           << "======" << endl
           << graph.extraMetadata() << endl
           << "======" << endl;
      cout << "Expecting the following extra metadata:" << endl
           << "======" << endl
           << expected_metadata << endl
           << "======" << endl;
      assert(graph.extraMetadata() == expected_metadata);

      if (first_pass) {
        cout << "Removing the index directory" << endl;
        graph.removeDumpedIndex(path);
        cout << endl;
      }
    } while (first_pass);

    cout << "Unfreezing the graph." << endl;
    graph.unfreeze();
    assert(!graph.frozen());
    cout << endl;

    cout << "Clearing the graph." << endl;
    graph.clear();
    check_graph_properties(graph, 0, 0, 0, false);
    cout << endl;
    cout << "Ensure that extra metadata are empty." << endl;
    assert(graph.extraMetadata().empty());

    cout << "Dumping the (empty) graph to " << path << " index directory which exists." << endl;
    graph.dump(path, true);
    // Dumping the graph freezes it.
    assert(graph.frozen());
    cout << endl;

    cout << "Unfreezing the graph." << endl;
    graph.unfreeze();
    assert(!graph.frozen());
    cout << endl;

    cout << "Loading the graph from the dumped (empty) index" << endl;
    graph.load(path);
    check_graph_properties(graph, 0, 0, 0, true);
    cout << endl;
    cout << "Ensure that extra metadata are still empty." << endl;
    assert(graph.extraMetadata().empty());

    cout << "Removing the index directory" << endl;
    graph.removeDumpedIndex(path);
    cout << endl;

    cout << "Removing the temporary directory (" << top_path << ")" << endl;
    rmdir(top_path);
  } else {
    cout << "Unable to get a temporary directory." << endl
         << "This test is not complete but is not considered as failed." << endl;
  }
  cout << endl;

  cout << "That's All, Folk!!!" << endl;
  return 0;

}
