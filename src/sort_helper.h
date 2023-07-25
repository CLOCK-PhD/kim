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

#ifndef __SORT_HELPER_H__
#define __SORT_HELPER_H__

#include <vector>     // for std::vector
#include <algorithm>  // for std::sort
#include <numeric>    // for std::iota
#include <utility>    // for std::swap
#include <functional> // for std::less
#include <cassert>    // for assert


namespace kim {

  /**
   * A swap template wrapper around the standard swap template
   * function.
   *
   * The standard swap template function used an explicit reference
   * for its parameters but is incompatible with symbolic reference
   * type (like in vector<bool>::reference for example).
   *
   * As an example, look the swap wrapper for
   * std::vector<bool>::reference below.
   *
   * \param t1 The first parameter to swap.
   *
   * \param t2 The second parameter to swap.
   */
  template <typename reference>
  void swap(reference t1, reference t2) {
    std::swap(t1, t2);
  }


  /**
   * Specialization of the swap wrapper template fucntion for
   * std::vector<bool>::reference parameters.
   *
   * \param t1 The first bit to swap.
   *
   * \param t2 The second bit to swap.
   */
  template <>
  void swap(std::vector<bool>::reference t1, std::vector<bool>::reference t2);

  /**
   * Helper class to sort a container using another container order.
   * The core code comes from:
   * https://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
   *
   * Elements of type T must be comparable using the Compare functor
   * template parameter.
   *
   * The Container must be having operator[](size_t) and size() methods.
   */
  template <typename T,
            typename Container = std::vector<T>,
            typename Compare = std::less<T> >
  class SortHelper {

  private:

    /**
     * The container of reference elements used to sort.
     */
    const Container &_ref;

    /**
     * The permutation used to sort any container.
     */
    std::vector<size_t> _permutation;

    /**
     * The comparison functor.
     */
    const Compare _compare;

  public:

    /**
     * Constructs an helper using the given container as reference.
     *
     * \param container The container used as a reference.
     */
    SortHelper(const Container &container):
      _ref(container), _permutation(), _compare() {
      reset();
    }

    /**
     * Reinitialize the permutation used to sort any container.
     *
     * This must be called if the reference container changes after the
     * helper has been created.
     */
    void reset() {
      _permutation.resize(_ref.size());
      std::iota(_permutation.begin(), _permutation.end(), 0);
      std::sort(_permutation.begin(), _permutation.end(), *this);
    }

    /**
     * Get the permutation used ot sort any container with this helper.
     *
     * \return Returns the permutation used ot sort any container with
     * this helper.
     */
    const std::vector<size_t> &permutation() const {
      return _permutation;
    }

    /**
     * Get the reference container used ot sort any container with this
     * helper.
     *
     * \return Returns the reference container used ot sort any
     * container with this helper.
     */
    const Container &reference() const {
      return _ref;
    }

    /**
     * This helper can be used as a comparison functor between two
     * positions (based on the reference container).
     *
     * \param i The position of the first element to compare in the
     * reference container.
     *
     * \param j The position of the second element to compare in the
     * reference container.
     *
     * \return Returns the (indirect) comparison of reference element
     * at rank i against reference element at rank j.
     */
    bool operator()(size_t i, size_t j) {
      return _compare(_ref[i], _ref[j]);
    }

    /**
     * Sort the given container using the reference container order.
     *
     * The elements of type reference must be swappable.
     *
     * The OtherContainer must be a container having operator[](size_t)
     * and size() methods and must be of the same size than the
     * reference container.
     *
     * \param container The container to sort.
     */
    template <typename U, typename OtherContainer = std::vector<U>, typename reference = U&>
    void sort(OtherContainer &container) const {
      sort<U, OtherContainer, reference>(container, _permutation);
    }

    /**
     * Sort the given container using the given permutation.
     *
     * The elements of type reference must be swappable.
     *
     * The OtherContainer must be a container having
     * operator[](size_t) and size() methods and must be of the same
     * size than the permutation vector.
     *
     * \param container The container to sort.
     *
     * \param permutation The sort order permutation to use.
     */
    template <typename U, typename OtherContainer = std::vector<U>, typename reference = U&>
    static void sort(OtherContainer &container, const std::vector<size_t> &permutation) {
      assert(container.size() == permutation.size());
      std::vector<bool> placed(permutation.size(), false);
      for (size_t i = 0; i < container.size(); ++i) {
        if (!placed[i]) {
          size_t actual_pos = i;
          size_t final_pos = permutation[i];
          while (i != final_pos) {
            // Placing the element at the actual position to its final
            // position.
            kim::swap<reference>(container[actual_pos], container[final_pos]);
            // This element at this position must not move anymore.
            placed[final_pos] = true;
            // Thus the element at the actual position (which was at the final
            // position) needs to be processed too.
            actual_pos = final_pos;
            final_pos = permutation[final_pos];
          }
          placed[i] = true;
        }
      }
    }

  };

}

#endif
// Local Variables:
// mode:c++
// End:

