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

#ifndef __BOUNDED_SIZE_STRING_H__
#define __BOUNDED_SIZE_STRING_H__

#include <cstdlib>
#include <string>
#include <iostream>

namespace kim {

  /**
   * Handles bounded size strings.
   *
   * This class mimics the string one but using a bounded size buffer to
   * store the string content.
   *
   * The aim of this class is to reduce the size of the object.
   */
  class BoundedSizeString {

  private:

    /**
     * The maximal size of all string instances (it can't be modified while some instance exists).
     */
    static size_t _maximal_size;

    /**
     * The number of instances.
     */
    static size_t _nb_instances;

    /**
     * The associated c string.
     */
    char *_str;

    /**
     * Copy the given c_str to the current bounded size string
     *
     * \param c_str The C string to copy.
     */
    void _copy(const char *c_str);

  public:

    /**
     * Try to set the maximal size to the given one for bounded size
     * string.
     *
     * The maximal size can be modified if and only there is no instance
     * of bounded size string.
     *
     * \remark This static method must be invoked at least once before
     * trying to create bounded size strings since by default, the
     * maximal size is set to 0.
     *
     * \param maximal_size The maximal size for the bounded size strings
     * to be created.
     *
     * \return Returns true if the size has been updated and false
     * otherwise.
     */
    static bool setMaximalSize(size_t maximal_size);

    /**
     * Get the maximal size of bounded size strings.
     *
     * \return Returns the maximal size of the bounded size strings.
     */
    inline static size_t getMaximalSize() {
      return _maximal_size;
    }

    /**
     * Get the number of running instances of bounded size strings.
     *
     * \return Returns the number of running instances of bounded size strings.
     */
    inline static size_t getNbInstances() {
      return _nb_instances;
    }

    /**
     * Build a bounded size string instance from the given C string and
     * updates the counter of number of instances.
     *
     * \param c_str The C string to copy. If the string length is
     * greater than the _maximal_size, it is truncated.
     */
    BoundedSizeString(const char *c_str = NULL);

    /**
     * Build a bounded size string instance from the given string and
     * updates the counter of number of instances.
     *
     * \param s The string to copy. If the string length is greater than
     * the _maximal_size, it is truncated.
     */
    BoundedSizeString(const std::string &s);

    /**
     * Copy contructor to build a bounded size string. It updates the
     * counter of number of instances.
     *
     * \param s The bounded size string to copy.
     */
    BoundedSizeString(const BoundedSizeString &s);

    /**
     * Free all allocated memory and updates the counter of number of
     * instance.
     */
    ~BoundedSizeString();

    /**
     * Assignment operator.
     *
     * \param s The bounded size string to assign (right operand).
     *
     * \return Return this bounded size string after assignment.
     */
    BoundedSizeString &operator=(const BoundedSizeString &s);

    /**
     * Assignment operator.
     *
     * \remark If the length of the right operand is greater than the
     * _maximal_size, it is truncated.
     *
     * \param s The string to assign (right operand).
     *
     * \return Return this bounded size string after assignment.
     */
    BoundedSizeString &operator=(const std::string &s);

    /**
     * Assignment operator.
     *
     * \remark If the length of the right operand is greater than the
     * _maximal_size, it is truncated.
     *
     * \param s The C string to assign (right operand).
     *
     * \return Return this bounded size string after assignment.
     */
    BoundedSizeString &operator=(const char *s);

    /**
     * Lexicographic comparison.
     *
     * \param s The bounded size string to compare.
     *
     * \return Return 1 if current bounded size string is
     * lexicographically strictly greater than the given one, 0 if the
     * two bounded size string are equals and -1 otherwise.
     */
    int compare(const BoundedSizeString &s) const;

    /**
     * right to left Lexicographic comparison.
     *
     * \warning Ending null character are taken into account in the
     * comparison. They are considered as normal characters which simply
     * have the lowest possible value.
     *
     * \param s The bounded size string to compare.
     *
     * \return Return 1 if current bounded size string is
     * lexicographically strictly greater than the given one when
     * reading from right to left, 0 if the two bounded size string are
     * equals and -1 otherwise.
     */
    int reverse_compare(const BoundedSizeString &s) const;

    /**
     * Comparison operator for equality.
     *
     * \param s The bounded size string to compare.
     *
     * \return Return true if both bounded size string are equals.
     */
    inline bool operator==(const BoundedSizeString &s) const {
      return compare(s) == 0;
    }

    /**
     * Comparison operator for inequality.
     *
     * \param s The bounded size string to compare.
     *
     * \return Return true if both bounded size string are different.
     */
    inline bool operator!=(const BoundedSizeString &s) const {
      return compare(s) != 0;
    }

    /**
     * Comparison operator for precedence.
     *
     * \param s The bounded size string to compare.
     *
     * \return Return true if current bounded size string
     * lexicographically strictly precedes the given bounded size
     * string.
     */
    inline bool operator<(const BoundedSizeString &s) const {
      return compare(s) < 0;
    }

    /**
     * Comparison operator for succession.
     *
     * \param s The bounded size string to compare.
     *
     * \return Return true if current bounded size string
     * lexicographically strictly succeeds the given bounded size
     * string.
     */
    inline bool operator>(const BoundedSizeString &s) const {
      return compare(s) > 0;
    }

    /**
     * Comparison operator for precedence.
     *
     * \param s The bounded size string to compare.
     *
     * \return Return true if current bounded size string
     * lexicographically precedes the given bounded size
     * string.
     */
    inline bool operator<=(const BoundedSizeString &s) const {
      return compare(s) <= 0;
    }

    /**
     * Comparison operator for succession.
     *
     * \param s The bounded size string to compare.
     *
     * \return Return true if current bounded size string
     * lexicographically succeeds the given bounded size
     * string.
     */
    inline bool operator>=(const BoundedSizeString &s) const {
      return compare(s) >= 0;
    }

    /**
     * Get the C string associated to the current bounded size string.
     *
     * \return Returns the C string associated to the current bounded
     * size string.
     */
    inline const char *c_str() const {
      return _str ? _str : "";
    }

    /**
     * Return the character at the given position in read-only mode.
     *
     * \param i The position of the wanted character (starting from
     * 0).
     *
     * \return Returns the (read-only) character at the given
     * position. If the position is greater than the maximal size, a
     * segmentation fault might (should) occur.
     */
    inline const char &operator[](size_t i) const {
      return _str[i];
    }

    /**
     * Return the character at the given position in read-write mode.
     *
     * \param i The position of the wanted character (starting from
     * 0).
     *
     * \return Returns the (read-write) character at the given
     * position. If the position is greater than the maximal size, a
     * segmentation fault might (should) occur.
     */
    inline char &operator[](size_t i) {
      return _str[i];
    }

    /**
     * Get the current bounded size string length.
     *
     * \warning Contrarily to the string objects, where length is
     * returned in constant time, the length is computed in linear time.
     *
     * \return Returns the current bounded size string length.
     */
    size_t length() const;

    /**
     * Reset the current bounded size string to the empty string.
     */
    void clear();

    /**
     * Test if the current bounded size string is empty.
     *
     * \return Returns true if the current bounded size string is empty.
     */
    bool empty() const;

    /**
     * Exchange the content of the given bounded size string and the
     * current one.
     *
     * \param s The bounded size string with which content is exchange.
     */
    inline void swap(BoundedSizeString &s) {
      std::swap(_str, s._str);
    }

  };

  // Some function overloading

  /**
   * Get line from stream into bounded size string.
   *
   * Extracts characters from is and stores them into s until the
   * delimitation character delim is found or the maximal length if s is
   * reached.
   *
   * The extraction also stops if the end of file is reached in is or if
   * some other error occurs during the input operation.
   *
   * If the delimiter is found, it is extracted and discarded (i.e. it
   * is not stored and the next input operation will begin after it).
   *
   * Note that any content in s before the call is replaced by the newly
   * extracted sequence.
   *
   * \param is The input stream from which characters are extracted.
   *
   * \param s The bounded size string to fill with the stream content.
   *
   * \param delim The delimiter used to end the stream extraction.
   *
   * \return Returns the stream after characters' extraction (included
   * the delimiter).
   */
  std::istream &getline(std::istream &is, BoundedSizeString& s, char delim = '\n');

  /**
   * Get the next visible string from stream into bounded size string.
   *
   * Discard any blank and non visible characters until one visible
   * character is found. Then extracts all visible (and non blank)
   * characters from is and stores them into s until a blanck ot
   * unvisible character is found or the maximal length if s is reached.
   *
   * The extraction also stops if the end of file is reached in is or if
   * some other error occurs during the input operation.
   *
   * Note that any content in s before the call is replaced by the newly
   * extracted sequence.
   *
   * \param is The input stream from which characters are extracted.
   *
   * \param s The bounded size string to fill with the stream content.
   *
   * \return Returns the stream after characters' extraction.
   */
  std::istream &operator>>(std::istream &is, BoundedSizeString &s);

  inline std::ostream &operator<<(std::ostream &os, const BoundedSizeString &s) {
    os << s.c_str();
    return os;
  }

  inline void swap(BoundedSizeString &s1, BoundedSizeString &s2) {
    s1.swap(s2);
  }

}

#endif
// Local Variables:
// mode:c++
// End:
