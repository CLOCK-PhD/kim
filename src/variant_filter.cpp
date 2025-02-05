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

#include "variant_filter.h"

#include "config.h"
#include "kim_exception.h"

#include <cassert>
#include <cstring>
#include <regex>

using namespace std;

BEGIN_KIM_NAMESPACE

static bool _noCaseCmp(const string &s1, const string &s2) {
  return strcasecmp(s1.c_str(), s2.c_str()) == 0;
}


static bool _stringOperation(VariantFilter::StringOperator op, const string &s1, const string &s2) {
  bool res = false;
  switch (op) {
  case VariantFilter::StringOperator::EQUAL: res = _noCaseCmp(s1, s2); break;
  case VariantFilter::StringOperator::REGEX: {
    regex r(s2);
    res = regex_match(s1, r); break;
  }
  case VariantFilter::StringOperator::NOT_EQUAL: res = !_noCaseCmp(s1, s2); break;
  }
  return res;
}

template <typename T>
static bool _numericalOperation(VariantFilter::NumericalOperator op, const T v1, const T v2) {
  bool res = false;
  switch (op) {
  case VariantFilter::NumericalOperator::EQUAL: res = v1 == v2; break;
  case VariantFilter::NumericalOperator::NOT_EQUAL: res = v1 != v2; break;
  case VariantFilter::NumericalOperator::LESS_THAN: res = v1 < v2; break;
  case VariantFilter::NumericalOperator::LESS_OR_EQUAL: res = v1 <= v2; break;
  case VariantFilter::NumericalOperator::GREATER_THAN: res = v1 > v2; break;
  case VariantFilter::NumericalOperator::GREATER_OR_EQUAL: res = v1 >= v2; break;
  }
  return res;
}

static bool _booleanOperation(VariantFilter::BooleanOperator op, const bool v1, const bool v2) {
  bool res = false;
  switch (op) {
  case VariantFilter::BooleanOperator::EQUAL: res = v1 == v2; break;
  case VariantFilter::BooleanOperator::NOT_EQUAL: res = v1 != v2; break;
  }
  return res;
}

VariantFilter::VariantFilter(vcfpp::BcfRecord &variant):
  _variant(variant)
{}

bool VariantFilter::onChrom(StringOperator op, const string &value) const {
  return _stringOperation(op, _variant.CHROM(), value);
}

bool VariantFilter::onPos(NumericalOperator op, const size_t value) const {
  return _numericalOperation<size_t>(op, _variant.POS(), value);
}

bool VariantFilter::onID(StringOperator op, const string &value) const {
  return _stringOperation(op, _variant.ID(), value);
}

bool VariantFilter::onQuality(NumericalOperator op, const double value) const {
  return _numericalOperation<double>(op, _variant.QUAL(), value);
}

bool VariantFilter::onFilter(StringOperator op, const string &value) const {
  return _stringOperation(op, _variant.FILTER(), value);
}

double to_double(const std::string s) {
  char *endptr;
  errno = 0;
  double v = strtod(s.c_str(), &endptr);
  if (errno) {
    throw strerror(errno);
  }
  if (*endptr) {
    Exception e;
    e << "Unable to parse numerical value '" << s << "'";
    throw e;
  }
  return v;
}

bool to_bool(const std::string s) {
  bool is_true = _noCaseCmp(s, "true") || _noCaseCmp(s, "1");
  bool is_false = _noCaseCmp(s, "false") || _noCaseCmp(s, "0");
  if (!is_true && !is_false) {
    Exception e;
    e << "Unable to parse boolean value '" << s << "'";
    throw e;
  }
  return is_true;
}

bool VariantFilter::onInfo(const string &key, ComparisonOperator op, const string &value) const {
  bool res = false;
  double scalar_value = 0;
  try {
    if (_variant.getINFO(key, scalar_value)) {
      double v = to_double(value);
      res = _numericalOperation(NumericalOperator(op), scalar_value, v);
    } else {
      string string_value;
      if (_variant.getINFO(key, string_value)) {
        res = _stringOperation(StringOperator(op), string_value, value);
      } else {
        string r1_str = "\\b";
        string r2_str;
        r1_str += key;
        r2_str = r1_str;
        r1_str += "\\b";
        r2_str += "=";
        regex r1(r1_str);
        string infos = _variant.allINFO();
        bool found = regex_search(infos, r1);
        if (found) {
          regex r2(r2_str + "=");
          if (regex_search(infos, r2)) {
            Exception e;
            e << "Unable to filter on key '" << key << "' (probably an array value)";
            throw e;
          }
        }
        bool bool_value = to_bool(value);
        res = _booleanOperation(BooleanOperator(op), found, bool_value);
      }
    }
  } catch (const invalid_argument &) {
    string infos = _variant.allINFO();
    if (op == ComparisonOperator::REGEX) {
      res = _stringOperation(StringOperator(op), "", value);
    } else {
      // try to convert value to numerical value
      bool universal_operator = ((op == ComparisonOperator::EQUAL)
                                 || (op == ComparisonOperator::NOT_EQUAL));
      try {
        if (!universal_operator) {
          // Can't be a boolean operator
          throw;
        }
        bool v = to_bool(value);
        res = _booleanOperation(BooleanOperator(op), 0., v);
      } catch (...) {
        try {
          double v = to_double(value);
          res = _numericalOperation(NumericalOperator(op), 0., v);
        } catch (...) {
          // Conversion failure
          // Falling back to string comparison
          if (!universal_operator) {
            // But it should not be compared
            throw;
          }
          res = _stringOperation(StringOperator(op), "", value);
        }
      }
    }
  }
  return res;
}

bool VariantFilter::onPloidy(NumericalOperator op, const size_t value) const {
  return _numericalOperation<size_t>(op, _variant.ploidy(), value);
}

bool VariantFilter::onSNP(BooleanOperator op, const bool value) const {
  return _booleanOperation(op, _variant.isSNP(), value);
}

bool VariantFilter::onMultiAllelicSNP(BooleanOperator op, const bool value) const {
  return _booleanOperation(op, _variant.isMultiAllelicSNP(), value);
}

bool VariantFilter::onSV(BooleanOperator op, const bool value) const {
  return _booleanOperation(op, _variant.isSV(), value);
}

bool VariantFilter::onIndel(BooleanOperator op, const bool value) const {
  return _booleanOperation(op, _variant.isIndel(), value);
}

ostream &operator<<(ostream &os, VariantFilter::ExpressionOperator op) {
  switch (op) {
  case VariantFilter::ExpressionOperator::AND: os << " AND "; break;
  case VariantFilter::ExpressionOperator::OR: os << " OR "; break;
  case VariantFilter::ExpressionOperator::NOT: os << " NOT "; break;
  case VariantFilter::ExpressionOperator::BLOCK_OPEN: os << "("; break;
  case VariantFilter::ExpressionOperator::BLOCK_CLOSE: os << ")"; break;
  }
  return os;
}

ostream &operator<<(ostream &os, VariantFilter::FilterOperation op) {
  switch (op) {
  case VariantFilter::FilterOperation::ON_CHROM: os << "CHROM"; break;
  case VariantFilter::FilterOperation::ON_POS: os << "POS"; break;
  case VariantFilter::FilterOperation::ON_ID: os << "ID"; break;
  case VariantFilter::FilterOperation::ON_QUAL: os << "QUAL"; break;
  case VariantFilter::FilterOperation::ON_FILTER: os << "FILTER"; break;
  case VariantFilter::FilterOperation::ON_INFO: os << "INFO"; break;
  case VariantFilter::FilterOperation::ON_PLOIDY: os << "PLOIDY"; break;
  case VariantFilter::FilterOperation::ON_SNP: os << "SNP"; break;
  case VariantFilter::FilterOperation::ON_MULTI_ALLELIC_SNP: os << "MSNP"; break;
  case VariantFilter::FilterOperation::ON_SV: os << "SV"; break;
  case VariantFilter::FilterOperation::ON_INDEL: os << "INDEL"; break;
  }
  return os;
}

ostream &operator<<(ostream &os, VariantFilter::ComparisonOperator op) {
  switch (op) {
  case VariantFilter::ComparisonOperator::EQUAL: os << "="; break;
  case VariantFilter::ComparisonOperator::REGEX: os << "~"; break;
  case VariantFilter::ComparisonOperator::NOT_EQUAL: os << "!="; break;
  case VariantFilter::ComparisonOperator::LESS_THAN: os << "<"; break;
  case VariantFilter::ComparisonOperator::LESS_OR_EQUAL: os << "<="; break;
  case VariantFilter::ComparisonOperator::GREATER_THAN: os << ">"; break;
  case VariantFilter::ComparisonOperator::GREATER_OR_EQUAL: os << ">="; break;
  }
  return os;
}

ostream &operator<<(ostream &os, VariantFilter::StringOperator op) {
  return os << VariantFilter::ComparisonOperator(op);
}

ostream &operator<<(ostream &os, VariantFilter::NumericalOperator op) {
  return os << VariantFilter::ComparisonOperator(op);
}

ostream &operator<<(ostream &os, VariantFilter::BooleanOperator op) {
  return os << VariantFilter::ComparisonOperator(op);
}

END_KIM_NAMESPACE
