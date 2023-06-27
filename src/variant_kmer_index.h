/******************************************************************************
*                                                                             *
*  Copyright © 2023      -- IGH / LIRMM / CNRS / UM                           *
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

#ifndef __VARIANT_KMER_INDEX_H__
#define __VARIANT_KMER_INDEX_H__

#include <cstdlib>
#include <list>
#include <vector>
#include <unordered_map>
#include <string>

#include <variant_identification.h>

namespace kim {

  /**
   * The Variant k-mer index
   */
  class VariantKmerIndex {


  public:

    /**
     * The type of the variant identifier.
     */
    typedef std::string VariantID_type;

    /**
     * The association between some given k-mer and some variant
     */
    struct VariantKmerAssociation {

      /**
       * The variant ID.
       */
      VariantID_type rs_id;

      /**
       * The rank of the k-mer in the variant k-mer markers.
       */
      size_t kmer_rank;

      /**
       * True if the k-mer appears natively in the genome.
       */
      bool in_genome;

    };

    /**
     * The type of the (partial) index associated to some k-mer
     * suffix.
     */
    typedef std::unordered_multimap<std::string, VariantKmerAssociation> PartialIndex_type;

  private:

    size_t k, k1, k2;

    std::vector<PartialIndex_type> index;

    bool warn;

    size_t checkFilenameCorrectness(const std::string &filename);

    void parseFile(const std::string &filename, size_t prefix, bool is_first);

  public:

    /**
     * Load the index of k-mers associated to variants using the files
     * form the given directory.
     *
     * \param path Directory where index files are stored.
     *
     * \param warn Activate (default) or deactivate warnings while
     * loading index.
     */
    VariantKmerIndex(const char *path, bool warn = true);

    /**
     * Get the length of the k-mers.
     *
     * \return Returns the length of the k-mers.
     */
    size_t getKmerLength() const;

    /**
     * Get the list of variant identifications involving the given
     * k-mer.
     *
     * \param kmer A string representing the k-mer (k-mer lookup is
     * case insensitive).
     *
     * \return Returns the list of variant identifications involving
     * the given k-mer.
     */
    std::list<VariantKmerAssociation> search(const std::string &kmer) const;

  };

}

#endif
// Local Variables:
// mode:c++
// End:
