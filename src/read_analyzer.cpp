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

#include "read_analyzer.h"

#include "config.h"
#include "gauss.h"

#include <cassert>
#include <cmath>
// #include <iostream>

using namespace std;

BEGIN_KIM_NAMESPACE

ReadAnalyzer::ReadAnalyzer(double alpha, double threshold, bool weak_mode):
  _variant_kmer_rates(), _variant_scores(), _weak_mode(weak_mode), _completed(false),
  _quantile(Gauss::quantile(alpha)),
  _analyzis(&ReadAnalyzer::_firstPassAnalyzis),
  alpha(alpha), threshold(threshold)
{
  assert(alpha >= 0);
  assert(alpha <= 1);
  assert(threshold >= 0);
  assert(threshold <= 1);
}

bool ReadAnalyzer::_firstPassAnalyzis(const ReadAnalyzer::VariantKmerRatesIteratorWrapper &vr, ReadAnalyzer::VariantStatistics &stats) {
  assert(!_completed);
  // update mean and variance using the Welford's algorithm
  float delta = vr.rate - stats.mean;
  ++stats.count;
  stats.mean += delta / stats.count;
  stats.variance += (vr.rate - stats.mean) * delta;
  return true;
}

bool ReadAnalyzer::_secondPassAnalyzis(const ReadAnalyzer::VariantKmerRatesIteratorWrapper &vr, ReadAnalyzer::VariantStatistics &stats) {
  assert(_completed);
  assert(stats.count > 0);
  Gauss g(stats.mean, sqrt(stats.variance));
  vr.zscore = g.getZscore(vr.rate);
  // double p = 1 - g.getCDF(vr.rate);
  // cerr << vr.node.variant << " rate " << vr.rate << " zscore " << vr.zscore << " => p = " << p << endl;
  return g.test(vr.rate, alpha);
}

const ReadAnalyzer::VariantKmerRates &ReadAnalyzer::analyze() {

  VariantKmerRates::iterator it = _variant_kmer_rates.begin();

  while (it != _variant_kmer_rates.end()) {
    VariantKmerRatesIteratorWrapper vr(it);
    // cerr << "looking for k-mers of variant '" << vr.node.variant << endl;
    vr.rate /= vr.node.in_degree;
    if (vr.rate > 1.0) vr.rate = 1.0;

    if (vr.rate >= threshold) {
      if ((this->*_analyzis)(vr, _variant_scores[vr.node])) {
        ++it;
      } else {
        it = _variant_kmer_rates.erase(it);
      }
    } else {
      // cerr << "Removing variant '" << vr.node.variant << " since rate is " << vr.rate << endl;
      it = _variant_kmer_rates.erase(it);
    }
  }

  return _variant_kmer_rates;

}

void ReadAnalyzer::reset() {
  _variant_kmer_rates.clear();
}

#ifndef __UNUSED__
# define __UNUSED__(x)
#endif

void ReadAnalyzer::add(const VariantNodesIndex::VariantNode &node, uint16_t __UNUSED__(kmer_rank), bool in_ref) {
  // cerr << "Add kmer at pos " << kmer_rank
  //      << " of variant " << node.variant
  //      << (in_ref
  //          ? (_weak_mode
  //             ? " (keepd since in ref but weak mode"
  //             : " (skipped since in ref and in strict mode")
  //          : "") << endl;
  if (!_completed || (_variant_scores.find(node) != _variant_scores.end())) {
    _variant_kmer_rates[node].rate += _weak_mode || !in_ref;
  }
  // cerr << "Now mean (count) = " << _variant_kmer_rates[node].rate << endl;
}

void ReadAnalyzer::complete() {
  if (!_completed) {
    // cerr << "Ending computation of variant scores" << endl;
    VariantScores::iterator it = _variant_scores.begin();
    while (it != _variant_scores.end()) {
      VariantScoreIteratorWrapper vs(it);
      double error = INFINITY;
      if (vs.stats.variance <= 0) {
        // This case occurs when vs.stats.count == 1 or
        // vs.node.in_degree == 1
        vs.stats.variance = 0;
      } else {
        // This division is the finalization of the Welford's
        // algorithm for the variance computation.
        vs.stats.variance /= vs.stats.count;
        // The division of the variance (in the sqrt()) comes from the
        // limit central theorem.
        error = _quantile * sqrt(vs.stats.variance / vs.stats.count);
        if (error < 0) error = -error;
      }
      if ((vs.stats.mean - error) < threshold) {
        // cerr << "The variant " << vs.node.variant << " is removed since it has "
        //      << vs.stats.mean << "±" << error << " k-mer ratio (from " << vs.stats.count << " reads)"
        //      << endl;
        it = _variant_scores.erase(it);
      } else {
        // cerr << "The variant " << vs.node.variant << " is kept since it has "
        //      << vs.stats.mean << "±" << error << " k-mer ratio (from " << vs.stats.count << " reads and variance " << vs.stats.variance << ")"
        //      << endl;
        ++it;
      }
    }
    _analyzis = &ReadAnalyzer::_secondPassAnalyzis;
    _completed = true;
  }
}

const ReadAnalyzer::VariantScores &ReadAnalyzer::result() {
  complete();
  return _variant_scores;
}

END_KIM_NAMESPACE
