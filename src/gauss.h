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

#ifndef __GAUSS_H__
#define __GAUSS_H__

#include <cstddef>
#include <vector>

namespace kim {

  /**
   * Helper class to handle Gaussian (normal) law.
   */
  class Gauss {

  private:

    /**
     * Values of the reduced centered Gaussian law for some probility
     * range are stored in a dedicated table.
     */
    struct _Table {
      /**
       * The value of the first stored probability
       */
      double start;

      /**
       * The step between two consecutive probability values.
       */
      double step;

      /**
       * The value of the probalities.
       */
      std::vector<double> values;

      /**
       * Test whether a given probability is stored in this table.
       *
       * \param v The probability to check.
       *
       * \return Returns true if the given probability is in this
       * table.
       */
      bool contains(double v) const;

      /**
       * Get the value corresponding to the cumulative density
       * function of the given value.
       *
       * \param v The probability to check.
       *
       * \return Returns the CDF of the given vlaue if stored in this
       * table or NaN otherwise.
       */
      double lookup(double v) const;

    };

    /**
     * The table for alpha starting from 0.50 to 0.90 by step of 0.01.
     */
    static const _Table _table_percent;

    /**
     * The table for alpha starting from 0.900 to 0.990 by step of
     * 0.005.
     */
    static const _Table _table_half_percent;

    /**
     * The table for alpha starting from 0.99 to °.999 by step of
     * 0.001.
     */
    static const _Table _table_1E3;

    /**
     * The table for alpha starting from 0.999 to 0.9999 by step
     * 10e-4.
     */
    static const _Table _table_1E4;

    /**
     * The table for alpha starting from 0.9999 to 0.99999 by step
     * 10e-5.
     */
    static const _Table _table_1E5;

    /**
     * The table for alpha starting from 0.99999 to 0.999999 by step
     * 10e-6.
     */
    static const _Table _table_1E6;

    /**
     * The table for alpha starting from 0.9999990 to 0.9999999 by
     * step 10e-7.
     */
    static const _Table _table_1E7;

    /**
     * The table for alpha starting from 0.9999999 to 0.99999999 by
     * step 10e-8.
     */
    static const _Table _table_1E8;

    /**
     * The table for alpha starting from 0.99999999 to 1 by step
     * 10e-9.
     */
    static const _Table _table_1E9;

    /**
     * The lowest quantile value for which the probability is assumed
     * to be equal ot 1.
     */
    static const double _upper_quantile;

    /**
     * The mean of this normal low handler.
     */
    const double _mu;

    /**
     * The standard deviation of this normal low handler.
     */
    const double _sigma;

  public:

    /**
     * Builds a Gaussian law handler.
     *
     * \param mu The mean of this normal law handler (default to 0).
     *
     * \param sigma The mean of this normal law handler (default to
     * 1).
     */
    inline Gauss(double mu = 0, double sigma = 1): _mu(mu), _sigma(sigma) {
    }

    /**
     * Get the Z-score associated to the given value.
     *
     * \param x The value for which to compute the Z-score.
     *
     * \return Returns the Z-score (the probability) of the given
     * value under this Gaussian law.
     */
    inline double getZscore(double x) const {
      return Zscore(x, _mu, _sigma);
    }

    /**
     * Get the Z-score associated to the given value following a
     * Gaussian law characterized by the given mean and standard
     * deviation.
     *
     * \param x The value for which to compute the Z-score.
     *
     * \param mu The mean of the Gaussian law.
     *
     * \param sigma The standard deviation of the Gaussian law.
     *
     * \return Returns the Z-score (the probability) of the given
     * value under the given Gaussian law parameters.
     */
    inline static double Zscore(double x, double mu, double sigma) {
      return (x - mu) / sigma;
    }

    /**
     * Get the quantile of the given probability.
     *
     * \param alpha The probability for which to compute the quantile.
     *
     * \return Returns the quantile of the cumulative distribution
     * frequency of the given probability under this Gaussian law.
     */
    inline double getQuantile(double alpha) const {
      return quantile(alpha, _mu, _sigma);
    }

    /**
     * Get the quantile of the given probability following a Gaussian
     * law characterized by the given mean and standard deviation.
     *
     * \param alpha The probability for which to compute the quantile.
     *
     * \param mu The mean of the Gaussian law.
     *
     * \param sigma The standard deviation of the Gaussian law.
     *
     * \return Returns the quantile the given probability under this
     * Gaussian law.
     */
    static double quantile(double alpha, double mu = 0, double sigma = 1);

    /**
     * Get the cumulative distribution frequency of the given
     * probability following this Gaussian law.
     *
     * \param x The probability for which to compute the CDF.
     *
     * \return Returns the cumulative distribution frequency of the
     * given probability under this Gaussian law. This is
     * \f$\Phi(-\infty, x)\f$ of \f$\mathcal{N}\left(\mu,
     * \sigma^2\right)\f$.
     */
    inline double getCDF(double x) const {
      return CDF(x, _mu, _sigma);
    }

    /**
     * Get the cumulative distribution frequency of the given
     * probability following the Gaussian law characterized by the
     * given mean and standard deviation.
     *
     * \param x The probability for which to compute the CDF.
     *
     * \param mu The mean of the Gaussian law (default to 0).
     *
     * \param sigma The standard deviation of the Gaussian law
     * (default to 1).
     *
     * \return Returns the cumulative distribution frequency of the
     * given probability under the Gaussian law characterized by the
     * given mean and standard deviation. This is \f$\Phi(-\infty,
     * x)\f$ of \f$\mathcal{N}\left(\mu, \sigma^2\right)\f$.
     */
    static double CDF(double x, double mu = 0, double sigma = 1);

  };

}

#endif
// Local Variables:
// mode:c++
// End:
