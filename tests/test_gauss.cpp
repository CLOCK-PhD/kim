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

#include <gauss.h>

#include <iostream>
#ifdef NDEBUG
#  undef NDEBUG
#endif
#include <cassert>

using namespace std;
using namespace kim;


void compare(double v1, double v2, double epsilon) {
  double d = v1 > v2 ? v1 - v2 : v2 - v1;
  cout << "|" << v1 << " - " << v2 << "| = " << d << " (and epsilon is " << epsilon << ")" << endl;
  assert(d <= epsilon);
}

int main() {

  for (int i = 0; i <= 10; ++i) {
    double x = (lrand48() % 80 - 40) / 10.;
    double p = Gauss::CDF(x);
    double x2 = Gauss::quantile(p);
    cout << "InvCDF(" << p << ") = " << x2 << " (expecting " << x << ")" << endl;
    compare(x, x2, 1e-1);
  }

  double values[][4] = {
    { -1.64,   0.05     , 1e-2},
    {  1.64,   0.95     , 1e-2},
    {  1.96,   0.975    , 1e-3},
    {  2.5758, 0.995    , 1e-4},
    {  3.0902, 0.999    , 1e-4},
    {  3.2905, 0.9995   , 1e-4},
    {  3.7190, 0.9999   , 1e-4},
    {  3.8906, 0.99995  , 1e-4},
    {  4.2649, 0.99999  , 1e-4},
    {  4.4172, 0.999995 , 1e-4}
  };

  for (auto val: values) {
    double q = val[0], p = val[1], e = val[2];
    double v = Gauss::CDF(q);
    cout << "CDF(" << q << ") = " << v << " (expecting " << p << ")" << endl;
    compare(v, p, e);
    v = Gauss::quantile(p);
    cout << "InvCDF(" << p << ") = " << v << " (expecting " << q << ")" << endl;
    compare(v, q, e);
  }

  double alpha_values[6] = { 0., 0.01, 0.025, 0.05, 0.1, 0.5 };
  double x_values[4] = {1., 1.64, 1.96, 2.57 };
  bool expected_results[3][12][8] = {
    /* LEFT_TAILED_HYPOTHESIS_TEST */
    {
      /* 0. */ {
        /* 1. */ true,
        /* -1. */ true,
        /* 1.64 */ true,
        /* -1.64 */ true,
        /* 1.96 */ true,
        /* -1.96 */ true,
        /* 2.57 */ true,
        /* -2.57 */ true
      },
      /* 1 - 0. */ {
        /* 1. */ false,
        /* -1. */ false,
        /* 1.64 */ false,
        /* -1.64 */ false,
        /* 1.96 */ false,
        /* -1.96 */ false,
        /* 2.57 */ false,
        /* -2.57 */ false
      },
      /* 0.01 */ {
        /* 1. */ true,
        /* -1. */ true,
        /* 1.64 */ true,
        /* -1.64 */ true,
        /* 1.96 */ true,
        /* -1.96 */ true,
        /* 2.57 */ true,
        /* -2.57 */ false
      },
      /* 1 - 0.01 */ {
        /* 1. */ false,
        /* -1. */ false,
        /* 1.64 */ false,
        /* -1.64 */ false,
        /* 1.96 */ false,
        /* -1.96 */ false,
        /* 2.57 */ true,
        /* -2.57 */ false
      },
      /* 0.025, */ {
        /* 1. */ true,
        /* -1. */ true,
        /* 1.64 */ true,
        /* -1.64 */ true,
        /* 1.96 */ true,
        /* -1.96 */ false,
        /* 2.57 */ true,
        /* -2.57 */ false
      },
      /* 1 - 0.025, */ {
        /* 1. */ false,
        /* -1. */ false,
        /* 1.64 */ false,
        /* -1.64 */ false,
        /* 1.96 */ true,
        /* -1.96 */ false,
        /* 2.57 */ true,
        /* -2.57 */ false
      },
      /* 0.05 */ {
        /* 1. */ true,
        /* -1. */ true,
        /* 1.64 */ true,
        /* -1.64 */ true,
        /* 1.96 */ true,
        /* -1.96 */ false,
        /* 2.57 */ true,
        /* -2.57 */ false
      },
      /* 1 - 0.05 */ {
        /* 1. */ false,
        /* -1. */ false,
        /* 1.64 */ false,
        /* -1.64 */ false,
        /* 1.96 */ true,
        /* -1.96 */ false,
        /* 2.57 */ true,
        /* -2.57 */ false
      },
      /* 0.1 */ {
        /* 1. */ true,
        /* -1. */ true,
        /* 1.64 */ true,
        /* -1.64 */ false,
        /* 1.96 */ true,
        /* -1.96 */ false,
        /* 2.57 */ true,
        /* -2.57 */ false
      },
      /* 1 - 0.1 */ {
        /* 1. */ false,
        /* -1. */ false,
        /* 1.64 */ true,
        /* -1.64 */ false,
        /* 1.96 */ true,
        /* -1.96 */ false,
        /* 2.57 */ true,
        /* -2.57 */ false
      },
      /* 0.5 */ {
        /* 1. */ true,
        /* -1. */ false,
        /* 1.64 */ true,
        /* -1.64 */ false,
        /* 1.96 */ true,
        /* -1.96 */ false,
        /* 2.57 */ true,
        /* -2.57 */ false
      },
      /* 1 - 0.5 */ {
        /* 1. */ true,
        /* -1. */ false,
        /* 1.64 */ true,
        /* -1.64 */ false,
        /* 1.96 */ true,
        /* -1.96 */ false,
        /* 2.57 */ true,
        /* -2.57 */ false
      },
    },

    /* RIGHT_TAILED_HYPOTHESIS_TEST */
    {
      /* 0. */ {
        /* 1. */ true,
        /* -1. */ true,
        /* 1.64 */ true,
        /* -1.64 */ true,
        /* 1.96 */ true,
        /* -1.96 */ true,
        /* 2.57 */ true,
        /* -2.57 */ true
      },
      /* 1 - 0. */ {
        /* 1. */ false,
        /* -1. */ false,
        /* 1.64 */ false,
        /* -1.64 */ false,
        /* 1.96 */ false,
        /* -1.96 */ false,
        /* 2.57 */ false,
        /* -2.57 */ false
      },
      /* 0.01 */ {
        /* 1. */ true,
        /* -1. */ true,
        /* 1.64 */ true,
        /* -1.64 */ true,
        /* 1.96 */ true,
        /* -1.96 */ true,
        /* 2.57 */ false,
        /* -2.57 */ true
      },
      /* 1 - 0.01 */ {
        /* 1. */ false,
        /* -1. */ false,
        /* 1.64 */ false,
        /* -1.64 */ false,
        /* 1.96 */ false,
        /* -1.96 */ false,
        /* 2.57 */ false,
        /* -2.57 */ true
      },
      /* 0.025, */ {
        /* 1. */ true,
        /* -1. */ true,
        /* 1.64 */ true,
        /* -1.64 */ true,
        /* 1.96 */ false,
        /* -1.96 */ true,
        /* 2.57 */ false,
        /* -2.57 */ true
      },
      /* 1 - 0.025, */ {
        /* 1. */ false,
        /* -1. */ false,
        /* 1.64 */ false,
        /* -1.64 */ false,
        /* 1.96 */ false,
        /* -1.96 */ true,
        /* 2.57 */ false,
        /* -2.57 */ true
      },
      /* 0.05 */ {
        /* 1. */ true,
        /* -1. */ true,
        /* 1.64 */ true,
        /* -1.64 */ true,
        /* 1.96 */ false,
        /* -1.96 */ true,
        /* 2.57 */ false,
        /* -2.57 */ true
      },
      /* 1 - 0.05 */ {
        /* 1. */ false,
        /* -1. */ false,
        /* 1.64 */ false,
        /* -1.64 */ false,
        /* 1.96 */ false,
        /* -1.96 */ true,
        /* 2.57 */ false,
        /* -2.57 */ true
      },
      /* 0.1 */ {
        /* 1. */ true,
        /* -1. */ true,
        /* 1.64 */ false,
        /* -1.64 */ true,
        /* 1.96 */ false,
        /* -1.96 */ true,
        /* 2.57 */ false,
        /* -2.57 */ true
      },
      /* 1 - 0.1 */ {
        /* 1. */ false,
        /* -1. */ false,
        /* 1.64 */ false,
        /* -1.64 */ true,
        /* 1.96 */ false,
        /* -1.96 */ true,
        /* 2.57 */ false,
        /* -2.57 */ true
      },
      /* 0.5 */ {
        /* 1. */ false,
        /* -1. */ true,
        /* 1.64 */ false,
        /* -1.64 */ true,
        /* 1.96 */ false,
        /* -1.96 */ true,
        /* 2.57 */ false,
        /* -2.57 */ true
      },
      /* 1 - 0.5 */ {
        /* 1. */ false,
        /* -1. */ true,
        /* 1.64 */ false,
        /* -1.64 */ true,
        /* 1.96 */ false,
        /* -1.96 */ true,
        /* 2.57 */ false,
        /* -2.57 */ true
      },
    },

    /* TWO_TAILED_HYPOTHESIS_TEST */
    {
      /* 0. */ {
        /* 1. */ true,
        /* -1. */ true,
        /* 1.64 */ true,
        /* -1.64 */ true,
        /* 1.96 */ true,
        /* -1.96 */ true,
        /* 2.57 */ true,
        /* -2.57 */ true
      },
      /* 1 - 0. */ {
        /* 1. */ false,
        /* -1. */ false,
        /* 1.64 */ false,
        /* -1.64 */ false,
        /* 1.96 */ false,
        /* -1.96 */ false,
        /* 2.57 */ false,
        /* -2.57 */ false
      },
      /* 0.01 */ {
        /* 1. */ true,
        /* -1. */ true,
        /* 1.64 */ true,
        /* -1.64 */ true,
        /* 1.96 */ true,
        /* -1.96 */ true,
        /* 2.57 */ true,
        /* -2.57 */ true
      },
      /* 1 - 0.01 */ {
        /* 1. */ false,
        /* -1. */ false,
        /* 1.64 */ false,
        /* -1.64 */ false,
        /* 1.96 */ false,
        /* -1.96 */ false,
        /* 2.57 */ false,
        /* -2.57 */ false
      },
      /* 0.025, */ {
        /* 1. */ true,
        /* -1. */ true,
        /* 1.64 */ true,
        /* -1.64 */ true,
        /* 1.96 */ true,
        /* -1.96 */ true,
        /* 2.57 */ false,
        /* -2.57 */ false
      },
      /* 1 - 0.025, */ {
        /* 1. */ false,
        /* -1. */ false,
        /* 1.64 */ false,
        /* -1.64 */ false,
        /* 1.96 */ false,
        /* -1.96 */ false,
        /* 2.57 */ false,
        /* -2.57 */ false
      },
      /* 0.05 */ {
        /* 1. */ true,
        /* -1. */ true,
        /* 1.64 */ true,
        /* -1.64 */ true,
        /* 1.96 */ false,
        /* -1.96 */ false,
        /* 2.57 */ false,
        /* -2.57 */ false
      },
      /* 1 - 0.05 */ {
        /* 1. */ false,
        /* -1. */ false,
        /* 1.64 */ false,
        /* -1.64 */ false,
        /* 1.96 */ false,
        /* -1.96 */ false,
        /* 2.57 */ false,
        /* -2.57 */ false
      },
      /* 0.1 */ {
        /* 1. */ true,
        /* -1. */ true,
        /* 1.64 */ true,
        /* -1.64 */ true,
        /* 1.96 */ false,
        /* -1.96 */ false,
        /* 2.57 */ false,
        /* -2.57 */ false
      },
      /* 1 - 0.1 */ {
        /* 1. */ false,
        /* -1. */ false,
        /* 1.64 */ false,
        /* -1.64 */ false,
        /* 1.96 */ false,
        /* -1.96 */ false,
        /* 2.57 */ false,
        /* -2.57 */ false
      },
      /* 0.5 */ {
        /* 1. */ false,
        /* -1. */ false,
        /* 1.64 */ false,
        /* -1.64 */ false,
        /* 1.96 */ false,
        /* -1.96 */ false,
        /* 2.57 */ false,
        /* -2.57 */ false
      },
      /* 1 - 0.5 */ {
        /* 1. */ false,
        /* -1. */ false,
        /* 1.64 */ false,
        /* -1.64 */ false,
        /* 1.96 */ false,
        /* -1.96 */ false,
        /* 2.57 */ false,
        /* -2.57 */ false
      },
    },

  };

  for (Gauss::HypothesisTest test = Gauss::LEFT_TAILED_HYPOTHESIS_TEST;
       (int) test <= (int) Gauss::TWO_TAILED_HYPOTHESIS_TEST;
       test = (Gauss::HypothesisTest) ((int) test + 1)) {
    for (size_t alpha_cpt = 0; alpha_cpt < 6; ++alpha_cpt) {
      for (size_t alpha_complement = 0; alpha_complement < 2; ++alpha_complement) {
        double alpha = alpha_complement ? 1 - alpha_values[alpha_cpt] : alpha_values[alpha_cpt];
        for (size_t x_cpt = 0; x_cpt < 4; ++x_cpt) {
          for (size_t x_opposite = 0; x_opposite < 2; ++x_opposite) {
            double x = x_opposite ? -x_values[x_cpt] : x_values[x_cpt];
            bool t;
            bool expected_t = expected_results[(int) test][2 * alpha_cpt + alpha_complement][2 * x_cpt + x_opposite];
            t = Gauss::hypothesisTest(x, alpha, test);
            cerr << "p(X < " << x << ") = " << Gauss::CDF(x) << endl;

            cout << "hypothesisTest(x = " << x << ","
                 << " alpha = " << alpha << ","
                 << " test = " << ((test == Gauss::LEFT_TAILED_HYPOTHESIS_TEST)
                                   ? "LEFT_TAILED_HYPOTHESIS_TEST"
                                   : ((test == Gauss::RIGHT_TAILED_HYPOTHESIS_TEST)
                                      ? "RIGHT_TAILED_HYPOTHESIS_TEST"
                                      : "TWO_TAILED_HYPOTHESIS_TEST")) << ") => "
                 << (t ? "true" : "false")
                 << " (expecting: " << (expected_t ? "true" : "false") << ")"
                 << endl;
            assert(t == expected_t);
          }
        }
      }
    }
  }

  return 0;
}
