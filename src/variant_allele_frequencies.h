/******************************************************************************
*                                                                             *
*  Copyright © 2025      -- IGH / LIRMM / CNRS / UM                           *
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

#ifndef __VARIANT_ALLELE_FREQUENCIES_H__
#define __VARIANT_ALLELE_FREQUENCIES_H__

#include <map>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include <variant_kmer_enumerator.h>

namespace kim {

  /**
   * Handles allele frequencies of some given variant
   *
   * For a given variant, this class allows to retrieve the allel
   * frequency for each of its alternative allele and for each given
   * population tag.
   */
  class VariantAlleleFrequencies {

  public:

    /**
     * Data structure to store allele frequencies of some variant
     * allele.
     */
    struct VariantPopulationAlleleFrequency {
      /**
       * The variant name (associated to one alternative allele).
       */
      std::string variant;
      /**
       * The allele frequency of each population tag.
       */
      std::map<std::string, float> population_af;

      /**
       * Handler for allele frequencies of some given variant.
       *
       * \param v The variant name.
       *
       * \param population_af The population allele frequencies for
       * this variant (default is empty).
       */
      VariantPopulationAlleleFrequency(const std::string &v,
                                       const std::map<std::string, float> &population_af = std::map<std::string, float>());

    };

  private:

    /**
     * The collection of alternative alleles with their frequencies
     * for each given population tag.
     */
    std::vector<VariantPopulationAlleleFrequency> _vpaf_array;

  public:

    /**
     * Build an allele frequencies handler for the given population tags.
     *
     * \param tags The population tags.
     */
    VariantAlleleFrequencies(const std::set<std::string> &tags);

    /**
     * The population tags.
     */
    const std::set<std::string> &tags;

    /**
     * Computes the allele frequencies of some variant for the handled
     * population tags.
     *
     * \param vke The variant k-mer enumerator for which to compute
     * the allele frequencies.
     *
     * \return Return the allele frequencies of the given variant for
     * the handled population tags.
     */
    VariantAlleleFrequencies &compute(const VariantKmerEnumerator &vke);

    /**
     * The collection of alternative alleles with their frequencies
     * for each handled population tag.
     *
     * \return Returns an array containing the alleles frequencies for
     * each handled populations of each alternative allele.
     */
    inline const std::vector<VariantPopulationAlleleFrequency> &getCurrentVariantAlleleFrequencies() const {
      return _vpaf_array;
    }

    /**
     * Get the variant name of each alternate allele.
     *
     * \return Returns the sorted list of the alternative alleles of
     * the last handled variant.
     */
    std::vector<std::string> getAlternateAlleleVariants() const;

    /**
     * Get the allele frequencies of each alternative allele for all
     * population tags.
     *
     * \param alternate_allele The name of the alternate allele.
     *
     * \return Returns the allele frequency of the given allele for
     * each population tag. If the given allele name doesn't belong to
     * the last handled variant, then an exception is thrown with an
     * explicit message.
     */
    const VariantPopulationAlleleFrequency &getAlleleFrequency(const std::string &alternate_allele) const;

    /**
     * Print the CSV header for the handled population tags
     *
     * \param os The stream on which to print the CSV header line.
     */
    void printCSVHeader(std::ostream &os) const;

    /**
     * Print the CSV row for the last handled variant.
     *
     * \param os The stream on which to print the CSV last handled
     * variant row.
     */
    void toCSV(std::ostream &os) const;
  };

  /**
   * Operator to print a variant allele frequency object.
   *
   * \param os The stream on which to print the CSV row of the last
   * variant handled by the given VariantAlleleFrequencies instance.
   *
   * \param vpaf The variant allele frequencies helper instance to
   * print.
   *
   * \return Returns the updated output stream.
   */
  std::ostream &operator<<(std::ostream& os, const VariantAlleleFrequencies &vaf);

}

#endif
// Local Variables:
// mode:c++
// End:
