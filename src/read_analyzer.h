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

#ifndef __READ_ANALYZER_H__
#define __READ_ANALYZER_H__

#include <cstddef>
#include <list>
#include <utility>
#include <map>

#include <variant_nodes_index.h>

namespace kim {

  /**
   * Helper class to analyse k-mer associated to variants in reads.
   */
  class ReadAnalyzer {

  public:

    /**
     * Type to store k-mer rates of variants in some read.
     */
    typedef std::map<VariantNodesIndex::VariantNode, double> VariantKmerRates;

    /**
     * A simple wrapper over VariantKmerRates::iterator to rename its
     * attributes.
     */
    struct VariantKmerRatesIteratorWrapper {

      /**
       * The node attribute.
       */
      const VariantNodesIndex::VariantNode& node;

      /**
       * The rate attribute.
       */
      double& rate;

      /**
       * Builds an iterator wrapper.
       *
       * \param it The iterator to wrap.
       */
      inline VariantKmerRatesIteratorWrapper(VariantKmerRates::iterator &it):
        node(it->first), rate(it->second) {
      }

    };

    /**
     * A simple wrapper over VariantKmerRates::const_iterator to
     * rename its attributes.
     */
    struct VariantKmerRatesConstIteratorWrapper {

      /**
       * The node attribute.
       */
      const VariantNodesIndex::VariantNode& node;

      /**
       * The rate attribute.
       */
      const double& rate;

      /**
       * Builds a const_iterator wrapper.
       *
       * \param it The const_iterator to wrap.
       */
      inline VariantKmerRatesConstIteratorWrapper(VariantKmerRates::const_iterator &it):
        node(it->first), rate(it->second) {
      }

    };

    /**
     * Statistics of the variant among all analyzed reads.
     */
    struct VariantStatistics {

      /**
       * The k-mer rate average.
       */
      float mean;

      /**
       * The k-mer rate variance.
       */
      float variance;

      /**
       * The number of analyzed reads having at least one (weak)
       * k-mer.
       */
      size_t count;

    };

    /**
     * Type to store k-mer statistics of variants in analyzed reads.
     */
    typedef std::map<VariantNodesIndex::VariantNode, VariantStatistics> VariantScores;

    /**
     * A simple wrapper over VariantScores::iterator to rename its
     * attributes.
     */
    struct VariantScoreIteratorWrapper {

      /**
       * The node attribute.
       */
      const VariantNodesIndex::VariantNode& node;

      /**
       * The statistics attribute.
       */
      VariantStatistics& stats;

      /**
       * Builds an iterator wrapper.
       *
       * \param it The iterator to wrap.
       */
      inline VariantScoreIteratorWrapper(VariantScores::iterator &it):
        node(it->first), stats(it->second) {
      }

    };

    /**
     * A simple wrapper over VariantScores::const_iterator to rename
     * its attributes.
     */
    struct VariantScoreConstIteratorWrapper {

      /**
       * The node attribute.
       */
      const VariantNodesIndex::VariantNode& node;

      /**
       * The statistics attribute.
       */
      const VariantStatistics& stats;

      /**
       * Builds a const_iterator wrapper.
       *
       * \param it The const_iterator to wrap.
       */
      inline VariantScoreConstIteratorWrapper(VariantScores::const_iterator &it):
        node(it->first), stats(it->second) {
      }

    };

  private:

    /**
     * The k-mer rates of variants found in the current analyzed read.
     */
    VariantKmerRates _variant_kmer_rates;

    /**
     * The detected variants statistics.
     */
    VariantScores _variant_scores;

    /**
     * Whether this analyzer uses all k-mers (true) or only k-mers
     * that are not in the reference (false).
     */
    const bool _weak_mode;

    /**
     * The analysis status.
     *
     * Once completed, no more informations can be added to this
     * analyzer.
     */
    bool _completed;

    /**
     * The quantile corresponding to alpha under the Centered Reduced
     * Gaussian law.
     */
    const double _quantile;

  public:

    /**
     * Builds a read analyzer, using the given significance and/or the
     * given threshold.
     *
     * \param alpha The type I error (significance) of the variant
     * analysis (the probability to reject a variant that really is in
     * the analyzed data). The type I error must be in the range [0;
     * 1].
     *
     * \param threshold The k-mer rate threshold to consider a variant
     * being present in some read. The threshold must be in the range
     * [0; 1].
     *
     * \param weak_mode Whether this analyzer uses all k-mers (true,
     * default) or only k-mers that are not in the reference (false).
     */
    ReadAnalyzer(double alpha, double threshold, bool weak_mode = true);

    /**
     * The type I error (significance) of the variant analysis (the
     * probability to reject a variant that really is in the analyzed
     * data).
     */
    const double alpha;

    /**
     * The threshold used to consider a variant in some read (thus in
     * the whole set of analyzed reads).
     */
    const double threshold;

    /**
     * Analyzer mode.
     *
     * A weak mode means that this analyzer uses all k-mers to [try
     * to] detect variants whereas a strict mode means that this
     * analyzer uses only k-mers that are not in the reference.
     *
     * \see strict()
     *
     * \return Returns true if this analyzer uses a weak mode.
     */
    inline bool weak() const {
      return _weak_mode;
    }

    /**
     * Analyzer mode.
     *
     * A weak mode means that this analyzer uses all k-mers to [try
     * to] detect variants whereas a strict mode means that this
     * analyzer uses only k-mers that are not in the reference.
     *
     * \see weak()
     *
     * \return Returns true if this analyzer uses a strict mode.
     */
    inline bool strict() const {
      return !_weak_mode;
    }

    /**
     * Perform the analysis of all collected variants since the last
     * time reset() was called or since the beginning.
     *
     * This method should be called each time a read has been fully
     * processed.
     *
     * \remark This method can't be called once the analysis is
     * completed.
     *
     * \return Returns the k-mer rates computed for all variants.
     */
    const VariantKmerRates &analyze();

    /**
     * Reset all collected variant (and associated informations).
     *
     * This method should be called before starting a new read
     * analysis.
     *
     * \remark This method can't be called once the analysis is
     * completed.
     */
    void reset();

    /**
     * Associate to the variant node a new k-mer position.
     *
     * For efficiency consideration, no test about duplicated k-mer
     * rank is performed. If the k-mer naturally exists in the
     * reference sequence(s), then it is not considered as a strict
     * k-mer.
     *
     * \param node The variant node.
     *
     * \param kmer_rank The k-mer rank to add.
     *
     * \param in_ref Existence of the k-mer in the reference
     * sequence(s).
     *
     * \remark This method can't be called once the analysis is
     * completed.
     */
    void add(const VariantNodesIndex::VariantNode &node, uint16_t kmer_rank, bool in_ref);

    /**
     * Get the analysis status (completed or not)
     *
     * \return Returns the analysis status.
     */
    inline bool completed() const {
      return _completed;
    }

    /**
     * Return the completed result.
     *
     * If the analysis wasn't already completed, then completes it
     * (ends the statistics computations) prior to returning the
     * analysis result.
     *
     * \return Returns the result of the complete analysis.
     */
    const VariantScores &result();

  };

}

#endif
// Local Variables:
// mode:c++
// End:
