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

#ifndef __VARIANT_NODES_INDEX_H__
#define __VARIANT_NODES_INDEX_H__

#include <cstdint>
#include <map>
#include <string>

namespace kim {

    /**
     * The type of each "variant" node of the bipartite graph.
     */
  class VariantNodesIndex: private std::map<std::string, uint16_t> {

  private:

    friend class KmerVariantEdgesSubindex;

  public:

    /**
     * The type describing the content of a variant node.
     */
    struct VariantNode {

      /**
       * The variant identifier.
       */
      const key_type &variant;

      /**
       * The node incoming degree.
       */
      mapped_type in_degree;

      /**
       * Constructs a VariantNode "view".
       *
       * \remark The "view" means that both the variant label and the
       * incoming degree used to build the view must exists while
       * current view is alive.
       *
       * \param it The read-only iterator "pointing to" the variant
       * node.
       */
      inline VariantNode(const_iterator it):
        variant(it->first),
        in_degree(it->second)
      {}

      /**
       * The undefined node with an empty identifier and a 0 incoming
       * degree.
       */
      static const VariantNode undefined;

    private:

      /**
       * Constructs a default VariantNode "view".
       */
      VariantNode();

    };

    using std::map<std::string, uint16_t>::iterator;
    using std::map<std::string, uint16_t>::const_iterator;
    using std::map<std::string, uint16_t>::cbegin;
    using std::map<std::string, uint16_t>::cend;
    using std::map<std::string, uint16_t>::crbegin;
    using std::map<std::string, uint16_t>::crend;
    using std::map<std::string, uint16_t>::size;
    using std::map<std::string, uint16_t>::max_size;
    using std::map<std::string, uint16_t>::empty;

    /**
     * Get the read-only node associated to the given variant.
     *
     * \param variant The variant for which the node is queried.
     *
     * \return Returns a "view" of the node associated to the given
     * variant if found ot the VariantNode::undefined node otherwise.
     */
    inline VariantNode getVariantNode(const std::string &variant) const {
      const_iterator it;
#ifdef _OPENMP
#pragma omp critical
#endif
      it = find(variant);
      return ((it == cend())
              ? VariantNode::undefined
              : VariantNode(it));
    }

    /**
     * Random access operator.
     *
     * This is a shortcut for getVariantNode() method.
     *
     * \param variant The variant for which the node is queried.
     *
     * \return Returns a "view" of the node associated to the given
     * variant if found ot the VariantNode::undefined node otherwise.
     */
    inline VariantNode operator[](const std::string &variant) const {
      return getVariantNode(variant);
    }

    /**
     * Get the number of k-mers associated to the given variant (the
     * incoming degree of its associated node).
     *
     * \param variant The variant for which the number of associated
     * k-mer is queried.
     *
     * \return Returns the number of k-mers associated to the given
     * variant.
     */
    inline uint16_t getVariantCount(const std::string &variant) const {
      return getVariantNode(variant).in_degree;
    }

    /**
     * Add (if necessary) a new node associated to the given variant.
     *
     * If there is no variant node for the given variant, a new node
     * is added with its incoming degree set to zero, otherwise it acts like find.
     *
     * \param variant The variant to add.
     *
     * \param increment_in_degree Increment the incoming degree of the
     * variant node.
     *
     * \return Returns a read-write iterator to the node corresponding
     * to the given variant.
     */
    iterator addVariantNode(const std::string &variant, bool increment_in_degree = false);

    /**
     * Add (if necessary) a new node associated to the given variant.
     *
     * This is a shortcut for addVariantNode() without incrementing
     * the incoming degree.
     *
     * \param variant The variant to add.
     *
     * \return Returns this index.
     */
    inline VariantNodesIndex &operator+=(const std::string &variant) {
      addVariantNode(variant);
      return *this;
    }

  };

}

#endif
// Local Variables:
// mode:c++
// End:
