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

#ifndef __VARIANT_FILTER_DRIVER_H__
#define __VARIANT_FILTER_DRIVER_H__

#include <string>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <vcfpp.h>
#pragma GCC diagnostic pop

#include <kim_exception.h>
#include <variant_filter.h>
#include <variant_filter_parser.h>

// Give Flex the prototype of yylex
#ifdef YY_DECL
#  undef YY_DECL
#endif
#define YY_DECL kim::VariantFilterParser::symbol_type yylex(kim::VariantFilterDriver &driver)

// Declare it for the parser's sake.
YY_DECL;

namespace kim {

  class VariantFilterScanner;

  /**
   * Exception associated to variant filter parse error.
   */
  class VariantFilterDriverException: public Exception {

  public:

    /**
     * Create a parse error exception associated to some variant
     * filter expression.
     *
     * \param driver The variant filter expression scanning and
     * parsing driver throwing this.
     */
    VariantFilterDriverException(const VariantFilterDriver &driver);

  };

  /**
   * Conducting the whole scanning and parsing of some filters to a
   * given Variant.
   */
  class VariantFilterDriver: public VariantFilter {

  private:

    /**
     * The filter expression string
     */
    std::string _filter;

    /**
     * The last (successfully) applied filter result.
     *
     * If no filter was applied, then this member is true.
     */
    bool _result;

    /**
     * The parsing error message if any
     */
    std::string _error_message;

    friend VariantFilterDriverException;
    friend VariantFilterParser;

  public:

    /**
     * Handles filters for a specific variant.
     *
     * \param variant The variant to apply filter on.
     */
    VariantFilterDriver(vcfpp::BcfRecord &variant);

    /**
     * Apply the given filter string to the handled variant.
     *
     * \param filter The filter expression to apply to the handled variant.
     *
     * \return Returns true if the variant successfully pass the given
     * filter and false otherwise. The returned value is also stored
     * in the result member attriubte. If filter expression is
     * syntaxically invalid or requires too much memory to be handled
     * (this situation should never occur), then an exception is
     * thrown.
     */
    bool apply(const std::string &filter);

    /**
     * Whether to generate parser debug traces or not.
     */
    bool debug_parsing;

    /**
     * Whether to generate scanner debug traces.
     */
    bool debug_scanning;

    /**
     * The token's location used by the scanner.
     */
    location scanner_location;

  };

}

#endif
// Local Variables:
// mode:c++
// End:
