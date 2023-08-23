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

#ifndef __KMER_NODES_SUBINDEX_H__
#define __KMER_NODES_SUBINDEX_H__

#include <vector>
#include <iostream>
#include <utility>
#include <bm.h>

#include <bounded_size_string.h>
#include <kim_exception.h>

namespace kim {

  /**
   * k-mer node sub-index exception.
   */
  class KmerNodesSubindexException: public Exception {

  public:

    /**
     * Create an exception dedicated to k-mer node sub-index.
     *
     * \param msg The initial message of this exception.
     */
    inline KmerNodesSubindexException(const std::string &msg = ""): Exception(msg) {}

  };


  /**
   * This class describes k-mer node collection sharing the same prefix.
   *
   * \see See class VariantKmerIndex.
   */
  class KmerNodesSubindex {

  public:

    /**
     * The type of each "k-mer" node of the bipartite graph.
     */
    struct KmerNode {

      /**
       * The suffix associated to the k-mer node.
       */
      BoundedSizeString suffix;

      /**
       * The information about the k-mer being in some of the
       * reference sequences.
       */
      bool in_reference;

      /**
       * Create a new k-mer node instance.
       *
       * \param suffix The associated suffix.
       *
       * \param in_reference The information about the k-mer being in
       * some of the reference sequences.
       */
      inline KmerNode(const BoundedSizeString &suffix = "", bool in_reference = false):
        suffix(suffix), in_reference(in_reference) {}

      /**
       * Output this k-mer node using JSON format.
       *
       * \see See the JavaScript Object Notation web page at
       * https://json.org/ for more informations.
       *
       * \param os The output stream.
       */
      inline void toJson(std::ostream &os) const {
        os << "{suffix:" << suffix << ",in_reference:" << in_reference << "}";
      }

      /**
       * Output this k-mer node using YAML format.
       *
       * \see See the Yet Another Markup Language web page at
       * https://yaml.org/ for more informations.
       *
       * \param os The output stream.
       */
      inline void toYaml(std::ostream &os) const {
        os << suffix << ": " << in_reference;
      }

      /**
       * Output this k-mer node using DSV format (Delimiter Separated
       * Value).
       *
       * \param os The output stream.
       *
       * \param delim The delimiter to separate the values.
       */
      inline void toDsv(std::ostream &os, char delim) const {
        os << suffix << delim << in_reference;
      }

    };

  private:

    /**
     * For a given prefix, k-mers are stored in the lexicographic
     */
    std::vector<BoundedSizeString> _suffixes;

    /**
     * Information about the k-mer existence in the reference.
     */
    bm::bvector<> _in_reference;

    /**
     * The sorting status of this sub-index.
     *
     * It must be sorted in order to use this sub-index and/or add
     * edges.
     */
    bool _sorted;

    /**
     * This sub-index can't be modified if _frozen is set to true.
     * Sub-index can't be queried if _frozen is set to false.
     */
    bool _frozen;

  public:

    /**
     * Initialize the subindex of k-mer nodes sharing the same prefix
     * (not stored here).
     *
     * \param estimated_nb_kmers Estimation of the number of k-mers
     * (used to preallocate memory and thus to reduce time due to
     * dynamic memory reallocation).
     *
     * \param sorted Consider future nodes of this sub-index to be
     * added in sorted order. If set to false, then each added k-mer
     * suffix is appended at the end even if duplicated. If true,
     * while k-mer suffixes are added in lexicographic order, then
     * duplicates are not appended.
     */
    explicit KmerNodesSubindex(size_t estimated_nb_kmers = 0, bool sorted = true);

    /**
     * Get the frozen state of this sub-index.
     *
     * \see See freeze() and unfreeze() methods.
     *
     * \return Returns true if this sub-index is frozen and false
     * otherwise.
     */
    inline bool frozen() const {
      return _frozen;
    }

    /**
     * Freeze this sub-index in order to prevent any further
     * modification and to allow queries.
     *
     * \remark Freezing the sub-index implies calling the sort()
     * method.
     *
     * \see See unfreeze() and frozen() methods.
     */
    void freeze();

    /**
     * Unfreeze this sub-index in order to allow modifications.
     *
     * Result of queries on unfrozen sub-index may lead to
     * inconsistent results and may lead to throwing
     * KmerNodesSubindexException.
     *
     * \see See freeze() and frozen() methods.
     */
    void unfreeze();

    /**
     * Get information about the sub-index being empty.
     *
     * \return Returns true if and only if this sub-index contains no
     * k-mer node.
     */
    inline bool empty() const {
      return _suffixes.empty();
    }

    /**
     * Get the number of indexed k-mers in this subindex.
     *
     * \return Returns the number of indexed k-mers in this subindex.
     */
    inline size_t size() const {
      return _suffixes.size();
    }

    /**
     * Get the reserved size for storing k-mer nodes in this subindex
     * to prevent dynamic reallocation.
     *
     * \return Returns the allocated space for k-mer nodes in this
     * subindex.
     */
    inline size_t capacity() const {
      return _suffixes.capacity();
    }

    /**
     * Reserve enough room to store n nodes.
     *
     * If the actual capacity is greater or equal to n, doesn't do
     * anything.
     *
     * \param n The capacity to reserve.
     */
    inline void reserve(size_t n) {
      if (frozen()) {
        throw KmerNodesSubindexException("Sub-index is frozen and calling the KmerNodesSubindex::reserve() method is not allowed.");
      }
      _suffixes.reserve(n);
    }

    /**
     * Shrink the memory in order to fit the data but have no extra
     * reserved memory.
     */
    inline void shrink_to_fit() {
      if (frozen()) {
        throw KmerNodesSubindexException("Sub-index is frozen and calling the KmerNodesSubindex::shrink_to_fit() method is not allowed.");
      }
      _suffixes.shrink_to_fit();
    }

    /**
     * Get a view of the k-mer node at the given rank.
     *
     * \warning No verification about node existence is performed,
     * thus this method may lead to unexpected results or cause a
     * segmentation fault.
     *
     * \param rank The rank of the k-mer node.
     *
     * \return Returns a view of the k-mer node at the given rank.
     */
    inline KmerNode operator[](size_t rank) const {
      return { _suffixes[rank], _in_reference[rank] };
    }

    /**
     * Return the first k-mer node of this sub-index.
     *
     * \warning If the sub-index is empty, this method may lead to
     * unexpected results or cause a segmentation fault.
     *
     * \return Returns a view of the first k-mer node of this sub-index.
     */
    inline KmerNode front() const {
      return operator[](0);
    }

    /**
     * Return the last k-mer node of this sub-index.
     *
     * \warning If the sub-index is empty, this method may lead to
     * unexpected results or cause a segmentation fault.
     *
     * \return Returns a view of the last k-mer node of this sub-index.
     */
    inline KmerNode back() const {
      return operator[](size() - 1);
    }

    /**
     * Add the given k-mer node to the sub-index.
     *
     * \param node The k-mer node to be added to this sub-index.
     *
     * \return Returns true if the node was added to this sub-index
     * and false if it already exists (be aware that the returned
     * value may be a false positive if nodes are not sorted).
     */
    bool add(const KmerNode &node);

    /**
     * Shortcut to create the given k-mer node.
     *
     * \param node The k-mer node to be added to this sub-index.
     *
     * \return Returns this sub-index.
     */
    inline KmerNodesSubindex &operator+=(const KmerNode &node) {
      add(node);
      return *this;
    }

    /**
     * Defines the status of the k-mer associated to the node at the
     * given rank.
     *
     * \warning This method throws a KmerNodesSubindexException with an
     * explicit message if the rank is out of bound.
     *
     * \param rank The rank of the k-mer node for which the status
     * must be set.
     *
     * \param status A true value means that the k-mer associated to
     * the node at the given rank appears naturally in some of the
     * reference sequences.
     */
    void setInReferenceKmer(size_t rank, bool status = true);

    /**
     * Defines the status of the k-mer associated with the given
     * suffix.
     *
     * \warning This method throws a KmerNodesSubindexException with an
     * explicit message if the sub-index is not sorted or if the k-mer
     * suffix hasn't been found in this sub-index.
     *
     * This method corresponds to first calling the getKmerNodeRank()
     * method then calling the preceeding setInReferenceKmer() method.
     *
     * \param suffix The suffix of the k-mer for which the status must
     * be set.
     *
     * \param status A true value means that the k-mer appears
     * naturally in some of the reference sequences.
     */
    inline void setInReferenceKmer(const BoundedSizeString &suffix, bool status = true) {
      return setInReferenceKmer(getKmerNodeRank(suffix), status);
    }

    /**
     * Get the status of the k-mer associated to the node at the given
     * rank.
     *
     * \warning This method throws a KmerNodesSubindexException with an
     * explicit message if the sub-index is not sorted or if the rank
     * is out of bound.
     *
     * \param rank The rank of the k-mer node for which the status is
     * asked.
     *
     * \return Returns true if the k-mer associated to the node at the
     * given rank in this sub-index appears in some of the reference
     * sequences.
     */
    bool isInReferenceKmer(size_t rank) const;

    /**
     * Get the status of the k-mer identified by the given suffix.
     *
     * \warning This method throws a KmerNodesSubindexException with an
     * explicit message if the sub-index is not sorted or if the k-mer
     * suffix hasn't been found in this sub-index.
     *
     * This method corresponds to first calling the getKmerNodeRank()
     * method then calling the preceeding isInReferenceKmer() method.
     *
     * \param suffix The suffix of the k-mer for which the status must
     * be asked.
     *
     * \return Returns true if the k-mer identified by the given
     * suffix appears in some of the reference sequences.
     */
    inline bool isInReferenceKmer(const BoundedSizeString &suffix) const {
      return isInReferenceKmer(getKmerNodeRank(suffix));
    }

    /**
     * Get the rank of the node associated to the given k-mer suffix
     * in the sub-index.
     *
     * This requires the sub-index to be sorted in order to perform a
     * dichotomic lookup. If not sorted, a KmerNodesSubindexException
     * is thrown with an explicit message.
     *
     * \param suffix The k-mer suffix associated to the k-mer node to
     * find.
     *
     * \return Returns the index of the k-mer node associated to the
     * given suffix in the this sub-index. If the k-mer suffix is not
     * found, then -1 is returned.
     */
    size_t getKmerNodeRank(const BoundedSizeString &suffix) const;

    /**
     * Get the rank of the given k-mer node in the sub-index.
     *
     * This requires the sub-index to be sorted in order to perform a
     * dichotomic lookup. If not sorted, a KmerNodesSubindexException
     * is thrown with an explicit message.
     *
     * \param node The k-mer node to find.
     *
     * \return Returns the index of the k-mer node in the this
     * sub-index. If the k-mer suffix is not found, then -1 is
     * returned.
     */
    inline size_t getKmerNodeRank(const KmerNode &node) const {
      return getKmerNodeRank(node.suffix);
    }

    /**
     * Get information about the k-mer nodes sorting order in this
     * sub-index.
     *
     * \return A true value means that the k-mer nodes are sorted in
     * the lexicographic order of the associated suffixes.
     */
    inline bool sorted() const {
      return _sorted;
    }

    /**
     * Sort the k-mer nodes in this sub-index.
     *
     * The sorting operation may throw a KmerNodesSubindexException if
     * two duplicated nodes (according to their associated suffix)
     * haven't the same informations about being in some reference
     * sequence.
     *
     * If the k-mer nodes are already sorted, does nothing.
     *
     * \return Returns the permutation order used to sort the
     * sub-index.
     */
    std::vector<size_t> sort();

    /**
     * Remove duplicated k-mer nodes in this sub-index.
     *
     * The unique method may throw a KmerNodesSubindexException if two
     * duplicated nodes (according to their associated suffix) haven't
     * the same informations about being in some reference sequence.
     *
     * \return return a bit vector such that '1' bits are associated
     * to nodes that remains in the sub-index and '0' bits to nodes
     * that were removed.
     */
    bm::bvector<> unique();

    /**
     * Expand current sub-index according to the given bit vector.
     *
     * After expansion, this sub-index will have size the size of the
     * schema.
     *
     * \remark The number of true bits in the given schema must
     * corresponds to the number of k-mer suffixes in this sub-index.
     *
     * \param schema Bit vector describing expansion to operate. The
     * i^th true bit corresponds to the position in the expanded
     * sub-index of the i^th k-mer suffix in the initial
     * sub-index. Each false bit after the i^th true bit and before
     * the (i+1)^th true bit corresponds to a duplicate of the i^th
     * k-mer suffix in the expanded sub-index.
     */
    void expand(const bm::bvector<> &schema);

    /**
     * Merges the given sub-index with this one.
     *
     * Both the given sub-index and this one must be sorted prior to
     * the merge operation. Otherwise a KmerNodesSubindexException is
     * thrown.
     *
     * \param idx The sub-index to merge with the current one.
     *
     * \param pos1 This bit vector is cleared and filled with the
     * final position of already present nodes in this sub-index (the
     * i-th bit set to '1' is at the position of the i-th node before
     * merging).
     *
     * \param pos2 This bit vector is cleared and filled with the
     * position of nodes of the right hand size sub-index that are in
     * this sub-index (the i-th bit set to '1' is at the position of
     * the i-th node of the right hand side).
     *
     * \return The current updated sub-index is returned.
     */
    KmerNodesSubindex &merge(const KmerNodesSubindex &idx, bm::bvector<> &pos1, bm::bvector<> &pos2);

    /**
     * Shortcut to merge the given index.
     *
     * \see See the merge() method.
     *
     * \param idx The sub-index to merge with the current one.
     *
     * \return The current updated sub-index is returned.
     */
    inline KmerNodesSubindex &operator+=(const KmerNodesSubindex &idx) {
      bm::bvector<> pos1, pos2;
      return merge(idx, pos1, pos2);
    }

    /**
     * Clear the sub-index.
     */
    void clear();

  };

  /**
   * Print the given k-mer node sub-index using JSON format.
   *
   * \see See KmerNodesSubindex::KmerNode::toJson() method for more
   * details.
   *
   * \param os The output stream.
   *
   * \param node The k-mer node to print.
   *
   * \return The output stream is returned to allow operator chaining.
   */
  inline std::ostream &operator<<(std::ostream &os, const KmerNodesSubindex &idx) {
    os << "[";
    for (size_t i = 0; i < idx.size(); ++i) {
      if (i) os << ",";
      idx[i].toJson(os);
    }
    os << "]";
    return os;
  }

  /**
   * Print the given k-mer node using JSON format.
   *
   * \see See KmerNodesSubindex::KmerNode::toJson() method for more
   * details.
   *
   * \param os The output stream.
   *
   * \param node The k-mer node to print.
   *
   * \return The output stream is returned to allow operator chaining.
   */
  inline std::ostream &operator<<(std::ostream &os, const KmerNodesSubindex::KmerNode &node) {
    node.toJson(os);
    return os;
  }

}

#endif
// Local Variables:
// mode:c++
// End:
