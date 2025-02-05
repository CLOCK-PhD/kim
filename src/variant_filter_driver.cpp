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

#include "variant_filter_driver.h"
#include "variant_filter_scanner.h"

#include "config.h"

#include <cassert>

using namespace std;

extern int yy_flex_debug;

BEGIN_KIM_NAMESPACE

VariantFilterDriverException::VariantFilterDriverException(const VariantFilterDriver &driver):Exception()
{
  const string filter = driver._filter;
  const location &loc = driver.scanner_location;
  _msg += "Filter expression parse error: ";
  _msg += driver._error_message;
  _msg += ".\n";
  bool prepend_line_numbers = ((loc.begin.line != 1) || (loc.end.line != 1));
  size_t nb_digits = (prepend_line_numbers ? to_string(loc.end.line).size() : 0);
  size_t first_pos = 0;
  size_t n = filter.size();
  int nb_l = 0;
  while (first_pos < n) {
    ++nb_l;
    if ((nb_l >= loc.begin.line) and (nb_l <= loc.end.line)) {
      _msg += ">";
    } else {
      _msg += " ";
    }
    _msg += " ";
    if (prepend_line_numbers) {
      string num = to_string(nb_l);
      assert(nb_digits >= num.size());
      _msg += string(nb_digits - num.size(), ' ');
      _msg += num;
      _msg += ": ";
    }
    size_t last_pos = filter.find_first_of("\n", first_pos);
    if (last_pos == string::npos) {
      last_pos = n;
    }
    assert(first_pos <= last_pos);
    assert(last_pos <= n);
    _msg += filter.substr(first_pos, last_pos - first_pos);
    _msg += "\n";
    if ((nb_l >= loc.begin.line) and (nb_l <= loc.end.line)) {
      size_t padding_length = nb_digits + 2 + 2 * prepend_line_numbers;
      size_t mark_length = 0;
      if (loc.begin.line == loc.end.line) {
        assert(loc.begin.column > 0);
        padding_length += loc.begin.column - 1;
        assert(loc.begin.column <= loc.end.column);
        mark_length += loc.end.column - loc.begin.column;
      } else {
        if (nb_l == loc.begin.line) {
          padding_length += loc.begin.column - 1;
          assert(last_pos - first_pos + 1 >= size_t(loc.begin.column));
          mark_length += last_pos - first_pos - loc.begin.column + 1;
        } else {
          if (nb_l == loc.end.line) {
            assert(loc.end.column > 0);
            mark_length += loc.end.column - 1;
          } else {
            assert(last_pos >= first_pos);
            mark_length += last_pos - first_pos;
          }
        }
      }
      mark_length += !mark_length; // mark_length is at least 1 ;-)
      assert(mark_length > 0);
      _msg += string(padding_length, ' ');
      _msg += string(mark_length, '^');
      _msg += "\n";
    }
    first_pos = last_pos + 1;
  }
}

VariantFilterDriver::VariantFilterDriver(vcfpp::BcfRecord &variant):
  VariantFilter(variant),
  _filter(), _result(true), _error_message(),
  debug_parsing(false), debug_scanning(false), scanner_location()
{}

bool VariantFilterDriver::apply(const std::string &filter) {
  _filter = filter;
  _result = true;
  _error_message.clear();
  scanner_location.begin.line = scanner_location.end.line = 1;
  scanner_location.begin.column = scanner_location.end.column = 1;
  YY_BUFFER_STATE state = yy_scan_string(filter.c_str());
  yy_flex_debug = debug_scanning;
  VariantFilterParser parse(*this);
  parse.set_debug_level(debug_parsing);
  int res = parse();
  yy_delete_buffer(state);
  if (res) {
    throw VariantFilterDriverException(*this);
  }
  return _result;
}

END_KIM_NAMESPACE
