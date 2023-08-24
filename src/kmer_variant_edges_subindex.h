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

#ifndef __KMER_VARIANT_EDGES_SUBINDEX_H__
#define __KMER_VARIANT_EDGES_SUBINDEX_H__

#include <cstdlib>
#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <bm.h>

#include <kmer_nodes_subindex.h>
#include <variant_nodes_index.h>

namespace kim {

  /**
   * Edges from k-mer nodes to variants sub-index exception.
   */
  class KmerVariantEdgesSubindexException: public Exception {

  public:

    /**
     * Create an exception dedicated to edges (from k-mer nodes to
     * variants) sub-index.
     *
     * \param msg The initial message of this exception.
     */
    inline KmerVariantEdgesSubindexException(const std::string &msg = ""): Exception(msg) {}

  };


  /**
   * The type of associations between k-mers suffixes that share a
   * common prefix and their associated variants (one k-mer may be
   * associated to several variants).
   */
  class KmerVariantEdgesSubindex {

  private:

    /**
     * The type of the association of some given k-mer with some variant.
     */
    struct _KmerVariantEdge {

      /**
       * Iterator pointing to the given variant in the index.
       */
      VariantNodesIndex::iterator variant_node_iterator;

      /**
       * The rank of the associated k-mer for the pointed variant.
       */
      uint16_t rank;

    };

  public:

    /**
     * The type of edge view.
     */
    struct Edge {

      /**
       * The starting k-mer suffix node view.
       */
      KmerNodesSubindex::KmerNode kmer_node;

      /**
       * The rank of the associated k-mer for the linked variant.
       */
      uint16_t rank;

      /**
       * The ending variant node view.
       */
      VariantNodesIndex::VariantNode variant_node;

    };

    /**
     * The type of the association of some given k-mer with some variant.
     */
    struct KmerVariantAssociation {

      /**
       * Iterator pointing to the given variant in the index.
       */
      VariantNodesIndex::VariantNode variant_node;

      /**
       * The rank of the associated k-mer for the pointed variant.
       */
      uint16_t rank;

      /**
       * Constructs a KmerVariantAssociation view from a _KmerVariantEdge.
       *
       * \param edge The edge between a k-mer node and a variant to
       * view.
       */
      inline KmerVariantAssociation(const _KmerVariantEdge &edge):
        variant_node(edge.variant_node_iterator),
        rank(edge.rank)
      {}

    };

  private:

    /**
     * The compressed bit vector (with rank/select capacity) where
     * each bit corresponds to one edge.
     *
     * A true value means that the corresponding edge is the first
     * one of some k-mer in the order of appearance.
     *
     * The select1(p) operation on this vector gives the position of
     * the first edge of the k-mer suffix at position p in the k-mer
     * nodes collection.
     *
     * \see This data structure comes from the BitMagic library
     * available under the Apache License Version 2.0
     * (http://www.apache.org/licenses/LICENSE-2.0). For more
     * information please visit: http://bitmagic.io
     */
    bm::bvector<> _first_kmer_edges;

    /**
     * The rank/select helper index to improve these queries.
     *
     * \see This data structure comes from the BitMagic library
     * available under the Apache License Version 2.0
     * (http://www.apache.org/licenses/LICENSE-2.0). For more
     * information please visit: http://bitmagic.io
     */
    bm::bvector<>::rs_index_type _rank_select_helper;

    /**
     * Edges (ordered by k-mers suffixes) from some k-mer to some
     * variant.
     *
     * The edge at index p in this vector corresponds to the k-mer
     * suffix given at index rank1(p) in the first_kmer_edges
     * compressed bit vector.
     */
    std::vector<_KmerVariantEdge> _kmer_variant_edges;

    /**
     * The k-mer nodes sub-index.
     */
    KmerNodesSubindex &_kmer_nodes;

    /**
     * The variant nodes index
     */
    VariantNodesIndex &_variant_nodes;

    /**
     * The maximum number of edges starting from some k-mer node.
     */
    size_t _max_associations;

    /**
     * This sub-index can't be modified if _frozen is set to true.
     * Sub-index can't be queried if _frozen is set to false.
     */
    bool _frozen;

    /**
     * Computes the maximum number of edges starting from a k-mer node
     * in this sub-index.
     *
     * \see See the getMaxAssociations() method.
     */
    void _updateMaxAssociations();

    /**
     * Get the position of the first edge of the specified k-mer node.
     *
     * \param pos The rank of the k-mer node in the k-mer node
     * subindex.
     *
     * \return Returns the rank of the first edge from the given k-mer
     * node or -1 if there is no edge (or if the node rank is not
     * valid).
    */
    size_t _selectKmerFirstEdgePosition(size_t pos) const;

  public:

    /**
     * Initialize the subindex of edges from k-mer nodes to variant
     * nodes sharing the same prefix (not stored here).
     *
     * \param kmer_nodes The k-mer nodes subindex.
     *
     * \param variant_nodes The variant nodes index.
     *
     * \param estimated_nb_edges Estimation of the number of edges
     * (used to preallocate memory and thus to reduce time due to
     * dynamic memory reallocation).
     */
    KmerVariantEdgesSubindex(KmerNodesSubindex &kmer_nodes,
                             VariantNodesIndex &variant_nodes,
                             size_t estimated_nb_edges = 0);

    /**
     * Get the related subindex of k-mer nodes.
     *
     * \return Returns the related subindex of k-mer nodes.
     */
    inline KmerNodesSubindex &kmerNodesSubindex() {
      return _kmer_nodes;
    }

    /**
     * Get the related read-only subindex of k-mer nodes.
     *
     * \return Returns the related read-only subindex of k-mer nodes.
     */
    inline const KmerNodesSubindex &kmerNodesSubindex() const {
      return _kmer_nodes;
    }

    /**
     * Get the related index of variant nodes.
     *
     * \return Returns the related index of variant nodes.
     */
    inline VariantNodesIndex &variantNodesIndex() {
      return _variant_nodes;
    }

    /**
     * Get the related read-only index of variant nodes.
     *
     * \return Returns the related read-only index of variant nodes.
     */
    inline const VariantNodesIndex &variantNodesIndex() const {
      return _variant_nodes;
    }

    /**
     * Compact this sub-index and its associated k-mer nodes sub-index.
     *
     * \remark Compacting this sub-index will call the
     * KmerNodesSubindex::sort() and KmerNodesSubindex::unique()
     * methods on the associated k-mer nodes sub-index if and only if
     * k-mer nodes sub-index is not already sorted and has the same
     * size than this sub-index. If the k-mer nodes sub-index is not
     * sorted has not the same size as this sub-index, an exception is
     * thrown with an explicit message. Be aware that the k-mer nodes
     * sub-index should not be sorted/made unique outside this method
     * if some edges have already been added to this sub-index. and
     * once compacted, only the merge() operation can be used to add
     * new edges.
     *
     * \see See unfreeze() and frozen() methods.
     */
    void compact();

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
     * \remark Freezing this sub-index also freezes the associated k-mer nodes
     * sub-index.
     *
     * \see See KmerNodesSubindex::freeze() method.
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
     * edges.
     */
    inline bool empty() const {
      return _kmer_variant_edges.empty();
    }

    /**
     * Get the number of indexed edges in this subindex.
     *
     * \return Returns the number of indexed edges in this subindex.
     */
    inline size_t size() const {
      return _kmer_variant_edges.size();
    }

    /**
     * Get the reserved size for storing edges in this subindex to
     * prevent dynamic reallocation.
     *
     * \return Returns the allocated space for edges in this subindex.
     */
    inline size_t capacity() const {
      return _kmer_variant_edges.capacity();
    }

    /**
     * Reserve enough room to store n edges.
     *
     * If the actual capacity is greater or equal to n, doesn't do
     * anything.
     *
     * \param n The capacity to reserve.
     */
    void reserve(size_t n);

    /**
     * Shrink the memory in order to fit the data but have no extra
     * reserved memory.
     */
    void shrink_to_fit();

    /**
     * Get a view of the edge at the given rank.
     *
     * \warning No verification about edge existence is performed,
     * thus this method may lead to unexpected results or cause a
     * segmentation fault.
     *
     * \param rank The rank of the edge.
     *
     * \return Returns a view of the edge at the given rank.
     */
    Edge operator[](size_t rank) const;

    /**
     * Return the first edge of this sub-index.
     *
     * \warning If the sub-index is empty, this method may lead to
     * unexpected results or cause a segmentation fault.
     *
     * \return Returns a view of the first edge of this sub-index.
     */
    inline Edge front() const {
      return operator[](0);
    }

    /**
     * Return the last edge of this sub-index.
     *
     * \warning If the sub-index is empty, this method may lead to
     * unexpected results or cause a segmentation fault.
     *
     * \return Returns a view of the last edge of this sub-index.
     */
    inline Edge back() const {
      return operator[](size() - 1);
    }

    /**
     * Get the list of variant associated to the given k-mer node.
     *
     * \param suffix The suffix associadted to he edges k-mer node source.
     *
     * \return Returns the list of all associations between the given
     * k-mer node and variant nodes.
     */
    std::list<KmerVariantAssociation> getKmerVariantAssociation(const BoundedSizeString &suffix) const;

    /**
     * Get the list of variant associated to the given k-mer node.
     *
     * \param kmer_node The edges k-mer node source.
     *
     * \return Returns the list of all associations between the given
     * k-mer node and variant nodes.
     */
    inline std::list<KmerVariantAssociation> getKmerVariantAssociation(KmerNodesSubindex::KmerNode kmer_node) const {
      return getKmerVariantAssociation(kmer_node.suffix);
    }

    /**
     * Add an edge between the given k-mer node to the given variant
     * node. The edge is labelled by the given rank.
     *
     * If the k-mer nodes sub-index is not sorted or if the given
     * k-mer node doesn't exist in the k-mer nodes sub-index, it is
     * added at the end of the k-mer nodes sub-index (possibly
     * breaking sub-index consistency). In such case, you should
     * consider using the merge method which will perform better.
     *
     * \see See the merge() method.
     *
     * \param node The k-mer node.
     *
     * \param variant The variant node.
     *
     * \param rank The k-mer rank within the k-mers associated to the
     * given variant. If the rank value is greater than the maximal
     * value (aka, 2^16 - 1), then a KmerVariantEdgesSubindexException
     * is thrown with an explicit message.
     *
     * \return Returns this sub-index.
     */
    KmerVariantEdgesSubindex &add(KmerNodesSubindex::KmerNode node,
                                  const std::string &variant,
                                  size_t rank);

    /**
     * Merges the given sub-index with this one.
     *
     * Both the given sub-index and this one must be compacted prior
     * to the merge operation. Otherwise a
     * KmerVariantEdgesSubindexException is thrown.
     *
     * \remark Both edges sub-indexes must share the same variant
     * index.
     *
     * \param idx The sub-index to merge with the current one.
     *
     * \return The current updated sub-index is returned.
     */
    KmerVariantEdgesSubindex &merge(KmerVariantEdgesSubindex &idx);

    /**
     * Shortcut to merge the given index.
     *
     * \see See the merge() method.
     *
     * \param idx The sub-index to merge with the current one.
     *
     * \return The current updated sub-index is returned.
     */
    inline KmerVariantEdgesSubindex &operator+=(KmerVariantEdgesSubindex &idx) {
      return merge(idx);
    }

    /**
     * Clear the sub-index.
     *
     * \remark This clears the k-mer nodes associated to this
     * sub-index, the edges (of course) of this sub-index and variant
     * nodes involved in edges that have now a 0 incoming degree.
     */
    void clear();

    /**
     * Get the maximal number of variant associations for a k-mer in
     * this sub-index.
     *
     * \return The maximum number of edges starting from some k-mer
     * node.
     */
    inline size_t getMaxAssociations() const {
      return _max_associations;
    }

    /**
     * Print the given edges (thus the involved nodes) to the given
     * output stream.
     *
     * \param os Output stream in which to print this sub-index.
     *
     * \param header If not empty, print the given string (by
     * prepending a '#' on each line of this header) at the beginning
     * of the given stream.
     *
     * \param footer If not empty, print the given string (by
     * prepending a '#' on each line of this footer) at the end of the
     * given stream.
     *
     * \param infos_in_header If true, print a summary in header
     * comment.
     *
     * \param legend_in_header If true, print column significations in
     * header comment.
     */
    void toStream(std::ostream &os,
                  const std::string &header = "Start of sub-index",
                  const std::string &footer = "End of sub-index",
                  bool infos_in_header = true,
                  bool legend_in_header = true) const;

  };

  /**
   * Print the given edges (thus the involved nodes) to the given
   * output stream.
   *
   * \param os Output stream in which to print this sub-index.
   *
   * \param idx The sub-index to print (right hand side).
   *
   * \return The output stream is returned.
   */
  inline std::ostream &operator<<(std::ostream &os, const KmerVariantEdgesSubindex &idx) {
    idx.toStream(os, "", "", true, false);
    return os;
  }

}

#endif
// Local Variables:
// mode:c++
// End:
