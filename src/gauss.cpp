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

#include "gauss.h"

#include "config.h"

#include <cmath>

using namespace std;

BEGIN_KIM_NAMESPACE

/* table for alpha starting from 0.50 to 0.90 by step of 0.01. */
const Gauss::_Table Gauss::_table_percent = {
  0.50, 0.01, {
    0.00000000 /* 0.50 */,
    0.02506891 /* 0.51 */,
    0.05015358 /* 0.52 */,
    0.07526986 /* 0.53 */,
    0.10043372 /* 0.54 */,
    0.12566135 /* 0.55 */,
    0.15096922 /* 0.56 */,
    0.17637416 /* 0.57 */,
    0.20189348 /* 0.58 */,
    0.22754498 /* 0.59 */,
    0.25334710 /* 0.60 */,
    0.27931903 /* 0.61 */,
    0.30548079 /* 0.62 */,
    0.33185335 /* 0.63 */,
    0.35845879 /* 0.64 */,
    0.38532047 /* 0.65 */,
    0.41246313 /* 0.66 */,
    0.43991317 /* 0.67 */,
    0.46769880 /* 0.68 */,
    0.49585035 /* 0.69 */,
    0.52440051 /* 0.70 */,
    0.55338472 /* 0.71 */,
    0.58284151 /* 0.72 */,
    0.61281299 /* 0.73 */,
    0.64334541 /* 0.74 */,
    0.67448975 /* 0.75 */,
    0.70630256 /* 0.76 */,
    0.73884685 /* 0.77 */,
    0.77219321 /* 0.78 */,
    0.80642125 /* 0.79 */,
    0.84162123 /* 0.80 */,
    0.87789630 /* 0.81 */,
    0.91536509 /* 0.82 */,
    0.95416525 /* 0.83 */,
    0.99445788 /* 0.84 */,
    1.03643339 /* 0.85 */,
    1.08031934 /* 0.86 */,
    1.12639113 /* 0.87 */,
    1.17498679 /* 0.88 */,
    1.22652812 /* 0.89 */,
    1.28155157 /* 0.90 */
  }
};

/* table for alpha starting from 0.900 to 0.990 by step of 0.005. */
const Gauss::_Table Gauss::_table_half_percent = {
  0.900, 0.005, {
    1.28155157 /* 0.900 */,
    1.31057911 /* 0.905 */,
    1.34075503 /* 0.910 */,
    1.37220381 /* 0.915 */,
    1.40507156 /* 0.920 */,
    1.43953147 /* 0.925 */,
    1.47579103 /* 0.930 */,
    1.51410189 /* 0.935 */,
    1.55477359 /* 0.940 */,
    1.59819314 /* 0.945 */,
    1.64485363 /* 0.950 */,
    1.69539771 /* 0.955 */,
    1.75068607 /* 0.960 */,
    1.81191067 /* 0.965 */,
    1.88079361 /* 0.970 */,
    1.95996398 /* 0.975 */,
    2.05374891 /* 0.980 */,
    2.17009038 /* 0.985 */,
    2.32634787 /* 0.990 */,
  }
};

/* table for alpha starting from 0.99 to °.999 by step of 0.001. */
const Gauss::_Table Gauss::_table_1E3 = {
  0.99, 1E-3, {
    2.32634787 /* 0.990 */,
    2.36561813 /* 0.991 */,
    2.40891555 /* 0.992 */,
    2.45726339 /* 0.993 */,
    2.51214433 /* 0.994 */,
    2.57582930 /* 0.995 */,
    2.65206981 /* 0.996 */,
    2.74778139 /* 0.997 */,
    2.87816174 /* 0.998 */,
    3.09023231 /* 0.999 */
  }
};

/* table for alpha starting from 0.999 to 0.9999 by step 10e-4. */
const Gauss::_Table Gauss::_table_1E4 = {
  0.999, 1E-4, {
    3.09023231 /* 0.999 */,
    3.12138915 /* 0.9991 */,
    3.15590676 /* 0.9992 */,
    3.19465105 /* 0.9993 */,
    3.23888012 /* 0.9994 */,
    3.29052673 /* 0.9995 */,
    3.35279478 /* 0.9996 */,
    3.43161440 /* 0.9997 */,
    3.54008380 /* 0.9998 */,
    3.71901649 /* 0.9999 */
  }
};

/* table for alpha starting from 0.9999 to 0.99999 by step 10e-5. */
const Gauss::_Table Gauss::_table_1E5 = {
  0.9999, 1E-5, {
    3.71901649 /* 0.99990 */,
    3.74554859 /* 0.99991 */,
    3.77501194 /* 0.99992 */,
    3.80816826 /* 0.99993 */,
    3.84612614 /* 0.99994 */,
    3.89059189 /* 0.99995 */,
    3.94440008 /* 0.99996 */,
    4.01281081 /* 0.99997 */,
    4.10747965 /* 0.99998 */,
    4.26489079 /* 0.99999 */
  }
};

/* table for alpha starting from 0.99999 to 0.999999 by step 10e-6. */
const Gauss::_Table Gauss::_table_1E6 = {
  0.99999, 1E-6, {
    4.26489079 /* 0.999990 */,
    4.28835654 /* 0.999991 */,
    4.31445102 /* 0.999992 */,
    4.34386118 /* 0.999993 */,
    4.37758785 /* 0.999994 */,
    4.41717341 /* 0.999995 */,
    4.46518391 /* 0.999996 */,
    4.52638932 /* 0.999997 */,
    4.61138236 /* 0.999998 */,
    4.75342430 /* 0.999999 */
  }
};

/* table for alpha starting from 0.9999990 to 0.9999999 by step 10e-7. */
const Gauss::_Table Gauss::_table_1E7 = {
  0.999999, 1E-7, {
    4.75342430 /* 0.9999990 */,
    4.77467242 /* 0.9999991 */,
    4.79832251 /* 0.9999992 */,
    4.82500464 /* 0.9999993 */,
    4.85563758 /* 0.9999994 */,
    4.89163846 /* 0.9999995 */,
    4.93536743 /* 0.9999996 */,
    4.99121712 /* 0.9999997 */,
    5.06895772 /* 0.9999998 */,
    5.19933752 /* 0.9999999 */
  }
};

/* table for alpha starting from 0.9999999 to 0.99999999 by step 10e-8. */
const Gauss::_Table Gauss::_table_1E8 = {
  0.9999999, 1E-8, {
    5.19933752 /* 0.99999990 */,
    5.21888852 /* 0.99999991 */,
    5.24066372 /* 0.99999992 */,
    5.26524823 /* 0.99999993 */,
    5.29349570 /* 0.99999994 */,
    5.32672377 /* 0.99999995 */,
    5.36712849 /* 0.99999996 */,
    5.41880097 /* 0.99999997 */,
    5.49085145 /* 0.99999998 */,
    5.61200065 /* 0.99999999 */
  }
};

/* table for alpha starting from 0.99999999 to 1 by step 10e-9. */
const Gauss::_Table Gauss::_table_1E9 = {
  0.99999999, 1E-9, {
    5.61200065 /* 0.999999990 */,
    5.63020007 /* 0.999999991 */,
    5.65047968 /* 0.999999992 */,
    5.67338789 /* 0.999999993 */,
    5.69972514 /* 0.999999994 */,
    5.73072757 /* 0.999999995 */,
    5.76845667 /* 0.999999996 */,
    5.81675557 /* 0.999999997 */,
    5.88419007 /* 0.999999998 */,
    5.99780047 /* 0.999999999 */,
    7.46780047 /* 1.000000000 */
  }
};

bool Gauss::_Table::contains(double v) const {
  // cerr << setprecision(20) << "Check if " << v << " is in table from " << start << " to " << (start + (step * (size - 1))) << endl;
  return (v >= start) && (v <= (start + (values.size() - 1) * step));
}

double Gauss::_Table::lookup(double v) const {
  if (!contains(v)) return NAN;
  size_t a = (v - start + (step / 2.)) / step;
  // cerr << v << " is in this table of size " << values.size() << " at rank " << a << endl;
  return values[a];
}

double Gauss::quantile(double alpha, double mu, double sigma) {
  double z = Zscore(alpha, mu, sigma);
  if (z < 0.5) return -quantile(1 - z);
  if (z > 1) return NAN;
  double res = _table_percent.lookup(z);
  if (!isnan(res)) return res;
  res = _table_half_percent.lookup(z);
  if (!isnan(res)) return res;
  res = _table_1E3.lookup(z);
  if (!isnan(res)) return res;
  res = _table_1E4.lookup(z);
  if (!isnan(res)) return res;
  res = _table_1E5.lookup(z);
  if (!isnan(res)) return res;
  res = _table_1E6.lookup(z);
  if (!isnan(res)) return res;
  res = _table_1E7.lookup(z);
  if (!isnan(res)) return res;
  res = _table_1E8.lookup(z);
  if (!isnan(res)) return res;
  res = _table_1E9.lookup(z);
  if (!isnan(res)) return res;
  return 1;
}

static const double SQRT2 = 1.4142135623730950488016887242096980785697L;

double Gauss::CDF(double x, double mu, double sigma) {
  double z = Zscore(x, mu, sigma);
  return erfc(-z / SQRT2) / 2;
}

END_KIM_NAMESPACE
