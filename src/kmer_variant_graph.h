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

#ifndef __KMER_VARIANT_GRAPH_H__
#define __KMER_VARIANT_GRAPH_H__

#include <cstdlib>
#include <filesystem>
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
    typedef std::vector<KmerVariantEdgesSubindex> KmerVariantEdgesIndex;

    /**
     * The type of Edge view.
     */
    struct Edge {

      /**
       * The source k-mer node label
       */
      const std::string &kmer;

      /**
       * The k-mer rank in the variant (edge label)
       */
      size_t rank;

      /**
       * The destination variant node label.
       */
      const std::string &variant;

    };

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
     * This graph can't be modified if _frozen is set to true. The
     * graph can't be queried if _frozen is set to false.
     */
    bool _frozen;

    /**
     * The extra metadata associated to this graph.
     */
    std::string _extra_metadata;

    /**
     * Clear and resize sub-indexes to the given size.
     *
     * \param total The number of sub-indexes.
     *
     * \param estimated_nb_kmers The number of estimated k-mers (to
     * preallocate enough room for each sub-index).
     *
     * \param sorted Set the node sub-indexes mode to this value.
     */
    void _resizeSubindexes(size_t total, size_t estimated_nb_kmers, bool sorted);

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
     * \param prefix The k-mer encoded prefix (the rank of the
     * sub-indexes to load).
     */
    void _parseFile(const std::filesystem::path &filename, size_t prefix);

    /**
     * parse the given metadata file in order to extract the content
     * of the 'Additional informations' field.
     */
    void _parseMetadata(const std::filesystem::path &filename);


  public:

    /**
     * Initialize the index of k-mers associated to variants using the
     * given parameters.
     *
     * \param settings The k-mer identification metric program
     * settings. The settings are frozen after this instance
     * construction and prior to any other initialization. If settings
     * index directory exist, then load the index from that
     * directory. In such case, the estimated_nb_kmers parameter is
     * unused.
     *
     * \param estimated_nb_kmers Estimation of the number of k-mers
     * (used to preallocate memory and thus to reduce time due to
     * dynamic memory reallocation). If not set, 1GB in the limit of
     * half the available memory is used. This parameter is used only
     * if the instance is not loaded from an existing dumped index.
     */
    KmerVariantGraph(Settings &settings, size_t estimated_nb_kmers = -1);

    /**
     * Load the index of k-mers associated to variants using the files
     * from the given directory.
     *
     * Any existing k-mer nodes, variant nodes and thus edges is
     * removed before loading the graph (by calling the clear() method
     * at the beginning).
     *
     * \remark Once loaded, this graph is frozen (see freeze(),
     * unfreeze() and frozen() methods).
     *
     * \param path Directory where index files to load are stored.
     *
     * \return Returns false if the index directory doesn't exist and
     * true if the directory exist and the index is successfully
     * loaded.
     */
    bool load(const std::filesystem::path &path);

    /**
     * Get this graph settings.
     *
     * \return Returns this graph settings.
     */
    inline const Settings &settings() const {
      return _settings;
    }

    /**
     * Get the number of indexed k-mers.
     *
     * \remark If graph is unfrozen, the result is not guaranteed and
     * might be overestimated.
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
     * Get the number of edges sub-indexes.
     *
     * \return Returns the number of sub-indexes.
     */
    inline size_t getNbSubindexes() const {
      return _edges.size();
    }

    /**
     * Get the variant nodes index.
     *
     * \return Returns the variant nodes index.
     */
    inline const VariantNodesIndex &getVariantNodesIndex() const {
      return _variant_nodes;
    }

    /**
     * Get the whole set of edges sub-indexes.
     *
     * \return Returns the whole set of edges sub-indexes.
     */
    inline const KmerVariantEdgesIndex &getKmerVariantEdgesIndex() const {
      return _edges;
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
     * Get the frozen state of this graph.
     *
     * \see See freeze() and unfreeze() methods.
     *
     * \return Returns true if this graph is frozen and false
     * otherwise.
     */
    inline bool frozen() const {
      return _frozen;
    }

    /**
     * Freeze this graph in order to prevent any further
     * modification and to allow queries.
     *
     * \see See unfreeze() and frozen() methods.
     */
    void freeze();

    /**
     * Unfreeze this graph in order to allow modifications.
     *
     * Result of queries on unfrozen graph may lead to inconsistent
     * results and may lead to throwing KmerVariantGraphException.
     *
     * \see See freeze() and frozen() methods.
     */
    void unfreeze();

    /**
     * Add an edge (thus its starting and ending nodes) from the given
     * k-mer to the given variant to this graph with the given rank.
     *
     * \param kmer The k-mer associated to the variant.
     *
     * \param rank The rank of the k-mer for this variant.
     *
     * \param variant The variant to associate to the given k-mer.
     *
     * \return Returns this graph.
     */
    KmerVariantGraph &add(const std::string &kmer, size_t rank, const  std::string &variant);

    /**
     * Another way to add an edge to the graph.
     *
     * \param e The edge to add.
     *
     * \return Returns this graph.
     */
    inline KmerVariantGraph &add(const Edge &e) {
      return add(e.kmer, e.rank, e.variant);
    }

    /**
     * A shortcut to add an edge to the graph.
     *
     * \param e The edge to add.
     *
     * \return Returns this graph.
     */
    inline KmerVariantGraph &operator+=(const Edge &e) {
      return add(e);
    }

    /**
     * (Un)set the given k-mer as a reference k-mer.
     *
     * This needs the graph to just have been frozen then unfrozen in
     * order to perform correct lookups in the k-mer sub-index.
     *
     * \see See the unit test test_kmer_variant_graph.cpp for an
     * example.
     *
     * \param kmer The k-mer to (un)set as reference.
     *
     * \param state The reference flag value to set.
     *
     * \return Returns true if the k-mer has been found (and (un)set
     * as expected).
     */
    bool setInReferenceKmer(const std::string &kmer, bool state = true);

    /**
     * Get the reference flag of the given k-mer.
     *
     * \param kmer The query k-mer.
     *
     * \return Returns true if the k-mer was found and is set as a
     * reference k-mer and false otherwise.
     */
    bool isInReferenceKmer(const std::string &kmer) const;

    /**
     * Get the number of k-mers associated to the given variant (the
     * incoming degree of its associated node).
     *
     * \param variant The variant for which the number of associated
     * k-mer is queried.
     *
     * \see See the VariantNodesIndex::getVariantCount() method.
     *
     * \return Returns the number of k-mers associated to the given
     * variant.
     */
    inline uint16_t getVariantCount(const std::string &variant) const {
      return _variant_nodes.getVariantCount(variant);
    }

    /**
     * Dump this index to the given directory.
     *
     * \remark The graph is frozen at the beginning of the dumping
     * operation.
     *
     * If the destination directory already exists and overwrite is
     * false, an exception with an explicit message is thrown.
     *
     * If the destination directory already exists and overwrite is
     * true, the existing dumped index should have all the files
     * associated to the k-mer nodes prefixes and no other ones and
     * may lead to unexpected behavior.
     *
     * \param path Directory where index files are stored.
     *
     * \param overwrite If true and if path directory already exists,
     * then invokes the removeDumpedIndex() static method.
     */
    void dump(const std::filesystem::path &path, bool overwrite = false);

    /**
     * Remove the dumped index located at the given directory and
     * corresponding to this graph settings (for the k-mer prefix
     * length only).
     *
     * \param path Directory where index files are stored.
     */
    void removeDumpedIndex(const std::filesystem::path &path) const;

    /**
     * Set extra metadata (used when index is dumped).
     *
     * \param metadata The metadata associated to the current graph.
     */
    void extraMetadata(const std::string &metadata);

    /**
     * Get extra metadata.
     *
     * When index is loaded from some directory, the value of extra
     * metadata is extracted from the 'Additional informations'
     * section (if it exists) of the metadata file. This value can be
     * set/overwritten by using the extraMetadata() method too.
     *
     * The metadata might be reformatted:
     * - Multiple paragraphs in informations might be reformatted in
     *   order to always start paragraphs by "- ".
     * - Multiple lines inside a parargraph might be reformatted in
     *   order to always start by two spaces.
     * - Empty lines will be removed.
     *
     * As an example, using the following string:
     *
     *     "Multiple paragraphs in informations might be reformatted in\n"
     *     "order to always start paragraphs by \"- \".\n\n"
     *     "Multiple lines inside a parargraph might be reformatted in\n"
     *     "order to always start by two spaces.\n\n\n\n"
     *     "- Empty lines will be removed."
     *
     * will be reformatted as if it was provided as:
     *
     *     "- Multiple paragraphs in informations might be reformatted in\n"
     *     "  order to always start paragraphs by \"- \".\n\n"
     *     "- Multiple lines inside a parargraph might be reformatted in\n"
     *     "  order to always start by two spaces.\n\n"
     *     "- Empty lines will be removed."
     *
     * \return Returns the metadata associated to the current graph.
     */
    inline const std::string &extraMetadata() const {
      return _extra_metadata;
    }

    /**
     * Clears this graph.
     */
    void clear();

  };


  /**
   * Get the given nucleotide as an integer.
   *
   * \param c The nucleotide to encode.
   *
   * \return Returns 0 for 'a' or 'A' (Adenine), 1 for 'c' or 'C'
   * (Cytosine), 2 for 'g' or 'G' (Guanine) and 3 for 't' or 'T'
   * (Thymine) or 'u' or 'U' (Uracil). For any other character, an
   * exception is thrown.
   */
  size_t encode(char c);

  /**
   * Encode the prefix of the given k-mer as an integer.
   *
   * \param kmer The k-mer for which the prefix is to be encoded.
   *
   * \param k1 The length of the prefix of the given k-mer. No
   * verification is done about the length of k1 which is supposed to
   * be less or equal than the length of the k-mer.
   *
   * \return Return the integer representation of the k-mer prefix of
   * length k1. Said differently, the returned value is the
   * lexicographic rank of the DNA string of length k1 among all DNA
   * string having size k1 (starting from 0).
   */
  size_t encode(const std::string &kmer, size_t k1);

  /**
   * Get the given integer as a nucleotide uppercase letter.
   *
   * \param v the encoded nucleotide.
   *
   * \return Returns 'A' if v is 0, 'C' if v is 1, 'G' if v is 2 and
   * 'T' if v is 3. For any other value of v, an exception is thrown.
   * (see encode(const char)) for more details.
   */
  char decode(size_t v);

  /**
   * Decode the given integer as a k1-mer.
   *
   * \param v The encoded value of the k1-mer to decode.  No
   * verification is done about the value of v which is supposed to be
   * less or equal than 4^k1.
   *
   * \param k1 The length of the k1-mer to decode.
   *
   * \return Return the k1-mer string corresponding to the given value
   * (see encode(const std::string &, size_t)) for more details.
   */
  std::string decode(size_t v, size_t k1);

  /**
   * Get a timestamp string corresponding to the current time.
   *
   * The defaut format requires at most 68 characters:
   * - 10 characters for the date (YYY-MM-DD)
   * - 1 characters for the 'T' letter between the date and the time
   * - 8 characters for the time (HH:MM:SS)
   * - 5 characters for the time zone (±0000)
   * - 23 fixed characters for the epoch time
   * - 20 characters max for the number of seconds (even if this
       program stay alive when the sun will become a red giant)
   * - 1 nul ending character
   * and thus, it does respect the length constraint.
   *
   * \param fmt The timestamp format (\see std::strftime() for format
   * specification). You must ensure that the given format will
   * produce a string that fits in a 256 character buffer.
   *
   * \return Returns a timestamp string corresponding to the current
   * time using the given format (default is the ISO 8601 format).
   */
  std::string timestamp(const char *fmt = "%FT%T%z (%s seconds since Epoch)");

}

#endif
// Local Variables:
// mode:c++
// End:
