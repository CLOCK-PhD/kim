/******************************************************************************
*                                                                             *
*  Copyright © 2024      -- IGH / LIRMM / CNRS / UM                           *
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

#ifndef __VARIANT_KMER_ENUMERATOR_H__
#define __VARIANT_KMER_ENUMERATOR_H__

#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <vcfpp.h>
#pragma GCC diagnostic pop

#include <dna_file_index.h>
#include <dna_file_reader.h>
#include <kim_settings.h>
#include <kmer_variant_graph.h>

namespace kim {

  /**
   * Enumerate k-mers associated to some variant.
   */
  class VariantKmerEnumerator {

  private:

    /**
     * Internal counter for missing variant ids.
     */
    static size_t _missing;

    /**
     * Name of the sequence of the last processed variant.
     */
    static std::string _chrom;

    /**
     * Name of the file havng the sequence of the last processed
     * variant.
     */
    static std::string _fname;

    /**
     * Internal reader to retrieve variant flanking sequences (of
     * length k-1).
     */
    static std::unique_ptr<DNAFileReader> _reader;

    /**
     * DNA file index used to retrieve efficiently the variant
     * flanking sequences.
     */
    static std::unique_ptr<DNAFileIndex> _index;

    /**
     * The program settings.
     */
    static std::unique_ptr<Settings> _settings;

    /**
     * The handled variant.
     */
    const vcfpp::BcfRecord &_v;

    /**
     * The variant ID. If the variant has no ID (ID is '.'), a new ID
     * is generated prefixed by 'MI' and suffixed by its unique number
     * (in ascending order of missing variants, starting from 1).
     */
    std::string _id;

    /**
     * The variant ID used when variant has multiple alleles. The ID
     * is built using the variant ID and by appending an underscore
     * ('_') then an allele counter.
     */
    std::string _sub_id;

    /**
     * The left flanking sequence of the variant. It is at most of
     * length (k - 1) but can be shorter if there is denegeracy
     * symbols or if the variant is located near the begin of the
     * sequence.
     */
    std::string _left_kmer;

    /**
     * The right flanking sequence of the variant. It is at most of
     * length (k - 1) but can be shorter if there is denegeracy
     * symbols or if the variant is located near the end of the
     * sequence.
     */
    std::string _right_kmer;

    /**
     * The current k-mer position (rank) in the "super-k-mer".
     */
    size_t _pos;

    /**
     * The variant alternate alleles
     */
    std::vector<std::string> _alt;

    /**
     * The position of the current allele in the _alt vector.
     */
    size_t _cur_alt;

    /**
     * The last computed k-mer for the current allele.
     */
    std::string _kmer;

  public:

    /**
     * Create a k-mer enumerator for the given variant.
     *
     * This uses the shared internal DNA file index that must be prior
     * computed using the static init() method.
     *
     * \param v The variant
     */
    VariantKmerEnumerator(const vcfpp::BcfRecord &v);

    /**
     * Compute the next k-mer associated to the variant handled by
     * this enumerator.
     *
     * This method is to be used in conjunction with the
     * getCurrentKmerVariantEdge() method. If some variant has more
     * than a unique alternate allele, then this method enumerates all
     * k-mers for all alleles in their order of appearance and adds a
     * suffix to the variant ID to distinguish them.
     *
     * \return Returns true is some k-mer was correctly computed and
     * false if there is no more available k-mer for the hnadled
     * variant.
     */
    bool nextVariantKmer();

    /**
     * Get the edge to add into the kmer-variant graph to link the
     * last computed k-mer with the handled variant.
     *
     * You must ensure that the last call to nextVariantKmer()
     * returned true in order to get the edge.
     *
     * \return Returns a kmer-variant graph edge from the last
     * computed k-mer to the variant handled by this enumerator. The
     * edge is labelled by the rank of the k-mer vor this variant.
     */
    KmerVariantGraph::Edge getCurrentKmerVariantEdge() const;

    /**
     * Get the DNA file index in read-only mode.
     *
     * \return Returns the reference to the read-only DNA file index.
     */
    inline static const DNAFileIndex &getIndex() {
      return *_index;
    }

    /**
     * Initialize the internal shared DNA file index needed to compute
     * the k-mers associated to some variant.
     *
     * This static method must be called prior to any use of
     * VariantKmerEnumerator instance.
     *
     * \param settings The program settings.
     *
     * \param dna_files The DNA files to index.
     *
     * \param bookmark_size The distance between bookmarks. Default is
     * every 1000bp, which allows to have an index of readonable size
     * without loosing too much time to retrieve variant location context.
     */
    static void init(const Settings &settings, const std::vector<std::string> &dna_files, size_t bookmark_size = 1000);

  };

}

#endif
// Local Variables:
// mode:c++
// End:
