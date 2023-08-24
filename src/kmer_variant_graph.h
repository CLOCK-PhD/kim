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

#ifndef __KMER_VARIANT_GRAPH_H__
#define __KMER_VARIANT_GRAPH_H__

#include <cstdlib>
#include <list>
#include <vector>
#include <string>

#include <bm.h>

#include <kim_settings.h>
#include <bounded_size_string.h>
#include <kmer_nodes_subindex.h>
#include <variant_nodes_index.h>
#include <kmer_variant_edges_subindex.h>

namespace kim {

  /**
   * The k-mer/Variant graph generic exception.
   */
  class KmerVariantGraphException: public Exception {

  public:

    /**
     * Create an exception dedicated to k-mer/variant graph.
     *
     * \param msg The initial message of this exception.
     */
    inline KmerVariantGraphException(const std::string &msg = ""): Exception(msg) {}

  };

  /**
   * The k-mer/Variant graph parse error.
   */
  class KmerVariantGraphParseError: public Exception {

  public:

    /**
     * Create an exception dedicated to k-mer/variant graph parse
     * error.
     *
     * \param msg The initial message of this exception.
     */
    inline KmerVariantGraphParseError(const std::string &msg = ""): Exception(msg) {}

  };

  /**
   * The k-mer/variant graph.
   *
   * This is mostly a bipartite graph where k-mers corresponds to the
   * first independent set of vertices and variants corresponds to the
   * second independent set of vertices. An edge in this graph
   * therefore corresponds to the association between a k-mer and a
   * variant.
   */
  class KmerVariantGraph {

  public:

    /**
     * The type of "k-mer" node collection.
     */
    typedef std::vector<KmerNodesSubindex> KmerNodeIndex;

    /**
     * The type of edges index.
     */
    typedef std::vector<KmerVariantEdgesSubindex>  KmerVariantEdgesIndex;

  private:

    /**
     * The index settings
     */
    Settings &_settings;

    /**
     * The total number of k-mer nodes in the graph.
     */
    size_t _nb_kmers;

    /**
     * The total number of edges in the graph.
     */
    size_t _nb_edges;

    /**
     * The array of all k-mer nodes sub-indexes.
     */
    KmerNodeIndex _kmer_nodes;

    /**
     * The set of all variant nodes.
     */
    VariantNodesIndex _variant_nodes;

    /**
     * The array of all k-mer to variant edges sub-indexes.
     */
    KmerVariantEdgesIndex _edges;

    /**
     * This graph can't be modified if _frozen is set to true.  Graph
     * can't be queried if _frozen is set to false.
     */
    bool _frozen;

    /**
     * Clear and resize sub-indexes to the given size.
     *
     * \param total The number of sub-indexes.
     *
     * \param estimated_nb_kmers The number of estimated k-mers (to
     * preallocate enough room for each sub-index).
     */
    void _resizeSubindexes(size_t total, size_t estimated_nb_kmers);

    /**
     * Ensure that the filename is correct (composed with exactly k_1
     * DNA symbols: A, C, G or T) and returns the endoded value for
     * this prefix file name.
     *
     * \param filename The name of the file containing an edges
     * sub-index (thus the k-mer nodes sub-index and some of the
     * variant nodes).
     *
     * \return Return the rank of the k-mer nodes and edges
     * sub-indexes.
     */
    size_t _checkFilenameCorrectness(const std::string &filename) const;

    /**
     * Parse the given file associated to the given k-mer encoded
     * prefix.
     *
     * If this file is the first analyzed one (i.e., settings are
     * unfrozen), then the prefix length, suffix length and then the
     * total length of the k-mers are computed from the index data
     * structure, then the settings are frozen. Otherwise, the
     * correctness of the index data structure is verified according
     * to the computed parameters.
     *
     * Said with more simplicity, only the first analyzed file must be
     * set this parameter to true and all other must be analyzed with
     * this parameter set to false.
     *
     * \param filename The full path of the file containing the given
     * sub-index to load.
     *
     * \param The k-mer encoded prefix (the rank of the sub-indexes to
     * load).
     */
    void _parseFile(const std::string &filename, size_t prefix);

  public:

    /**
     * Initialize the index of k-mers associated to variants using the
     * given directory and the given parameters.
     *
     * \param settings The k-mer identification metric program
     * settings. The settings are frozen after this instance
     * construction and prior to any other initialization.
     *
     * \param estimated_nb_kmers Estimation of the number of k-mers
     * (used to preallocate memory and thus to reduce time due to
     * dynamic memory reallocation). If not set, 1GB in the limit of
     * half the available memory is used.
     */
    KmerVariantGraph(Settings &settings, size_t estimated_nb_kmers = -1);

    /**
     * Load the index of k-mers associated to variants using the files
     * from the given directory.
     *
     * \param path Directory where index files are stored.
     *
     * \param settings The k-mer identification metric program
     * settings. Notice that only the warning status will be
     * unmodified. After instance construction, the settings will be
     * updated and frozen.
     */
    KmerVariantGraph(const char *path, Settings &settings);

    /**
     * Get the number of indexed k-mers.
     *
     * \return Returns the number of indexed k-mers.
     */
    inline size_t getNbKmers() const {
      return _nb_kmers;
    }

    /**
     * Get the number of indexed variants.
     *
     * \return Returns the number of indexed variants.
     */
    inline size_t getNbVariants() const {
      return _variant_nodes.size();
    }

    /**
     * Get the number of associations between k-mers and variants.
     *
     * \return Returns the number of associations between k-mers and
     * variants (the number of edges of the bipartite graph).
     */
    inline size_t getNbKmerVariantEdges() const {
      return _nb_edges;
    }


    /**
     * Get the list of variant involving the given k-mer.
     *
     * \param kmer A string representing the k-mer (k-mer lookup is
     * case insensitive).
     *
     * \return Returns the list of variant identifications involving
     * the given k-mer.
     */
    std::list<KmerVariantEdgesSubindex::KmerVariantAssociation> search(const std::string &kmer) const;


    /**
     *
     */
    bool addKmerNode(const std::string &kmer, bool in_reference);

    void freeze();


    /**
     *
     */
    KmerVariantGraph &add(const std::string &variant, const std::string &kmer, size_t rank);

    /**
     *
     */
    bool setInReferenceKmer(const std::string &kmer, bool state = true);

    bool isInReferenceKmer(const std::string &kmer) const;

    inline uint16_t getVariantCount(const std::string &variant) const {
      return _variant_nodes.getVariantCount(variant);
    }

    /**
     * Dump this index to the given directory.
     *
     * \param path Directory where index files are stored.
     */
    void dump(const char *path);

  };

}

#endif
// Local Variables:
// mode:c++
// End:
