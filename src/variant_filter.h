/******************************************************************************
*                                                                             *
*  Copyright © 2024-2025 -- IGH / LIRMM / CNRS / UM                           *
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

#ifndef __VARIANT_FILTER_H__
#define __VARIANT_FILTER_H__

#include <iostream>
#include <string>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <vcfpp.h>
#pragma GCC diagnostic pop

namespace kim {

  /**
   * Tests whether a given variant meets some conditions.
   */
  class VariantFilter {

  public:

    /**
     * Logical expression operators.
     */
    enum class ExpressionOperator {
      AND,         /**<< Conjunction */
      OR,          /**<< Disjunction */
      NOT,         /**<< Negation */
      BLOCK_OPEN,  /**< Opening priority block */
      BLOCK_CLOSE  /**< Closing priority block */
    };

    /**
     * Variant filter operation.
     */
    enum class FilterOperation {
      ON_CHROM,              /**< Operation based on 'CHROM' VCF field */
      ON_POS,                /**< Operation based on 'POS' VCF field */
      ON_ID,                 /**< Operation based on 'ID' VCF field */
      ON_QUAL,               /**< Operation based on 'QUAL' VCF field */
      ON_FILTER,             /**< Operation based on 'FILTER' VCF field */
      ON_INFO,               /**< Operation based on 'INFO' VCF field */
      ON_PLOIDY,             /**< Operation based on ploidy property */
      ON_SNP,                /**< Operation based on single nucleotide polymorphism property. */
      ON_MULTI_ALLELIC_SNP,  /**< Operation based on multiple single nucleotide polymorphism property. */
      ON_SV,                 /**< Operation based on structural variant property. */
      ON_INDEL               /**< Operation based on insertion/deletion property. */
    };

    /**
     * Generic comparison operator.
     */
    enum class ComparisonOperator {
      EQUAL,             /**< Test for equality operator */
      REGEX,             /**< Test for regular expression matching operator */
      NOT_EQUAL,         /**< Test for inequality operator */
      LESS_THAN,         /**< Test for strict precedence operator */
      LESS_OR_EQUAL,     /**< Test for previousness operator */
      GREATER_THAN,      /**< Test for strict nextness operator */
      GREATER_OR_EQUAL   /**< Test for nextness operator */
    };

    /**
     * Subset of generic comparison operator for string comparison.
     */
    enum class StringOperator {
      EQUAL = int(ComparisonOperator::EQUAL),         /**< Test for equality operator */
      REGEX = int(ComparisonOperator::REGEX),         /**< Test for regular expression matching operator */
      NOT_EQUAL = int(ComparisonOperator::NOT_EQUAL)  /**< Test for inequality operator */
    };

    /**
     * Subset of generic comparison operator for numerical comparison.
     */
    enum class NumericalOperator {
      EQUAL = int(ComparisonOperator::EQUAL),                       /**< Test for equality operator */
      NOT_EQUAL = int(ComparisonOperator::NOT_EQUAL),               /**< Test for equality operator */
      LESS_THAN = int(ComparisonOperator::LESS_THAN),               /**< Test for strict precedence operator */
      LESS_OR_EQUAL = int(ComparisonOperator::LESS_OR_EQUAL),       /**< Test for precedence operator */
      GREATER_THAN = int(ComparisonOperator::GREATER_THAN),         /**< Test for strict nextness operator */
      GREATER_OR_EQUAL = int(ComparisonOperator::GREATER_OR_EQUAL)  /**< Test for nextness operator */
    };

    /**
     * Subset of generic comparison operator for boolean comparison.
     */
    enum class BooleanOperator {
      EQUAL = int(ComparisonOperator::EQUAL),         /**< Test for equality operator */
      NOT_EQUAL = int(ComparisonOperator::NOT_EQUAL)  /**< Test for inequality operator */
    };

  private:

    /**
     * The variant to apply some filters.
     */
    vcfpp::BcfRecord &_variant;

  public:

    /**
     * Handles filters on the given variant.
     *
     * \param variant The variant to handle.
     */
    VariantFilter(vcfpp::BcfRecord &variant);

    /**
     * Filter based on the 'CHROM' VCF field.
     *
     * \par Example:
     *
     *     vf.onChrom(StringOperator::REGEX, "seq[1-3]");
     *
     * \param op The operator to apply between the variant 'CHROM'
     * field and the given value.
     *
     * \param value The value to compare with the variant 'CHROM'
     * field using the given comparison operator.
     *
     * \return Returns the result of the comparison.
     */
    bool onChrom(StringOperator op, const std::string &value) const;

    /**
     * Filter based on the 'POS' VCF field.
     *
     * \par Example:
     *
     *     vf.onPos(NumericalOperator::LESS_THAN, 12345);
     *
     * \param op The operator to apply between the variant 'POS' field
     * and the given value.
     *
     * \param value The value to compare with the variant 'POS' field
     * using the given comparison operator.
     *
     * \return Returns the result of the comparison.
     */
    bool onPos(NumericalOperator op, const size_t value) const;

    /**
     * Filter based on the 'ID' VCF field.
     *
     * \par Example:
     *
     *     vf.onID(StringOperator::EQUAL, "RS12345");
     *
     * \param op The operator to apply between the variant 'ID' field
     * and the given value.
     *
     * \param value The value to compare with the variant 'ID' field
     * using the given comparison operator.
     *
     * \return Returns the result of the comparison.
     */
    bool onID(StringOperator op, const std::string &value) const;

    /**
     * Filter based on the 'QUAL' VCF field.
     *
     * \par Example:
     *
     *     vf.onQuality(StringOperator::GREATER_OR_EQUAL, 23);
     *
     * \param op The operator to apply between the variant 'QUAL'
     * field and the given value.
     *
     * \param value The value to compare with the variant 'QUAL' field
     * using the given comparison operator.
     *
     * \return Returns the result of the comparison.
     */
    bool onQuality(NumericalOperator op, const double value) const;

    /**
     * Filter based on the 'FILTER' VCF field.
     *
     * \par Example:
     *
     *     vf.onFilter(StringOperator::EQUAL, "PASS");
     *
     * \param op The operator to apply between the variant 'FILTER'
     * field and the given value.
     *
     * \param value The value to compare with the variant 'FILTER'
     * field using the given comparison operator.
     *
     * \return Returns the result of the comparison.
     */
    bool onFilter(StringOperator op, const std::string &value) const;

    /**
     * Filter based on some given key of the 'INFO' VCF field.
     *
     * \par Example:
     *
     *     vf.onFilter("1000G", StringOperator::EQUAL, "true")
     *     || vf.onFilter("DP", StringOperator::GREATER_OR_EQUAL, "154");
     *     || vf.onFilter("VC", StringOperator::REGEX, "([mM]?)[sS][nN][vVpP]");
     *
     * \param key The key (from the INFO field) to comapre with the
     * given value using the given comparison operator.
     *
     * \param op The operator to apply between the value of the given
     * key in the variant 'INFO' field and the given value.
     *
     * \param value The value to compare with the value of the key
     * (from the variant 'INFO' field) using the given comparison
     * operator.
     *
     * \return Returns the result of the comparison.
     */
    bool onInfo(const std::string &key, ComparisonOperator op, const std::string &value) const;

    /**
     * Filter based on the variant ploidy.
     *
     * \par Example:
     *
     *     vf.onPloidy(NumericalOperator::EQUAL, 2)
     *
     * \param op The operator to apply between the variant ploidy and
     * the given value.
     *
     * \param value The value to compare with the variant ploidy using
     * the given comparison operator.
     *
     * \return Returns the result of the comparison.
     */
    bool onPloidy(NumericalOperator op, const size_t value) const;

    /**
     * Filter based on the variant SNP status.
     *
     * \par Example:
     *
     *     vf.onSNP(NumericalOperator::EQUAL, true)
     *
     * \remark A Multi allelic SNP is not considered as a SNP by this
     * filter (see onMultiAllelicSNP).
     *
     * \param op The operator to apply between the variant SNP status
     * and the given value.
     *
     * \param value The value to compare with the variant SNP status
     * using the given comparison operator.
     *
     * \return Returns the result of the comparison.
     */
    bool onSNP(BooleanOperator op, const bool value) const;

    /**
     * Filter based on the variant multi-allelic SNP status.
     *
     * \par Example:
     *
     *     vf.onMultiAllelicSNP(NumericalOperator::EQUAL, true)
     *
     * \param op The operator to apply between the variant
     * multi-allelic SNP status and the given value.
     *
     * \param value The value to compare with the variant
     * multi-allelic SNP status using the given comparison operator.
     *
     * \return Returns the result of the comparison.
     */
    bool onMultiAllelicSNP(BooleanOperator op, const bool value) const;

    /**
     * Filter based on the variant SV (structural variation) status.
     *
     * \par Example:
     *
     *     vf.onSV(NumericalOperator::NOT_EQUAL, false)
     *
     * \param op The operator to apply between the variant SV status
     * and the given value.
     *
     * \param value The value to compare with the variant SV status
     * using the given comparison operator.
     *
     * \return Returns the result of the comparison.
     */
    bool onSV(BooleanOperator op, const bool value) const;

    /**
     * Filter based on the variant Indel (instertion/deletion) status.
     *
     * \par Example:
     *
     *     vf.onSV(NumericalOperator::NOT_EQUAL, false)
     *
     * \param op The operator to apply between the variant Indel
     * status and the given value.
     *
     * \param value The value to compare with the variant Indel status
     * using the given comparison operator.
     *
     * \return Returns the result of the comparison.
     */
    bool onIndel(BooleanOperator op, const bool value) const;

  };

  /**
   * Pretty print expression operator.
   *
   * \param os The output stream on which to print the expression
   * operator.
   *
   * \param op The expression operator to print.
   *
   * \return Returns the updated stream.
   */
  std::ostream &operator<<(std::ostream &os, VariantFilter::ExpressionOperator op);

  /**
   * Pretty print filter operation.
   *
   * \param os The output stream on which to print the filter
   * operation.
   *
   * \param op The filter operation to print.
   *
   * \return Returns the updated stream.
   */
  std::ostream &operator<<(std::ostream &os, VariantFilter::FilterOperation op);

  /**
   * Pretty print generic comparison operator.
   *
   * \param os The output stream on which to print the comparison
   * operator.
   *
   * \param op The comparison operator to print.
   *
   * \return Returns the updated stream.
   */
  std::ostream &operator<<(std::ostream &os, VariantFilter::ComparisonOperator op);

  /**
   * Pretty print string comparison operator.
   *
   * \param os The output stream on which to print the comparison
   * operator.
   *
   * \param op The comparison operator to print.
   *
   * \return Returns the updated stream.
   */
  std::ostream &operator<<(std::ostream &os, VariantFilter::StringOperator op);

  /**
   * Pretty print numerical comparison operator.
   *
   * \param os The output stream on which to print the comparison
   * operator.
   *
   * \param op The comparison operator to print.
   *
   * \return Returns the updated stream.
   */
  std::ostream &operator<<(std::ostream &os, VariantFilter::NumericalOperator op);

  /**
   * Pretty print boolean comparison operator.
   *
   * \param os The output stream on which to print the comparison
   * operator.
   *
   * \param op The comparison operator to print.
   *
   * \return Returns the updated stream.
   */
  std::ostream &operator<<(std::ostream &os, VariantFilter::BooleanOperator op);

}

#endif
// Local Variables:
// mode:c++
// End:
