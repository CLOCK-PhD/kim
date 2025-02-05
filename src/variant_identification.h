/******************************************************************************
*                                                                             *
*  Copyright © 2023-2025 -- IGH / LIRMM / CNRS / UM                           *
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

#ifndef __VARIANT_IDENITIFICATION_H__
#define __VARIANT_IDENITIFICATION_H__

#include <string>
#include <map>
#include <list>

namespace kim {

  /**
   * Handles reads associated to some variation.
   */
  class VariantIdentification {

  public:

    /**
     * The type of the read identifier.
     */
    typedef std::string ReadID_type;

    /**
     * The type of the association between reads and k-mer (that
     * identify some variation) positions.
     */
    typedef std::map<VariantIdentification::ReadID_type, std::list<size_t> > VariantReadKmersAssoc;

  private:

    /**
     * The positions of k-mers associated to variants in reads.
     */
    VariantReadKmersAssoc _informations;

  public:

    /**
     * Builds a variant identification object with no informations.
     */
    VariantIdentification();

    /**
     * Add some information associated to some variant.
     *
     * This method consider that a new k-mer (at the given position)
     * is associated to some potential variant for the specified read
     * ID.
     *
     * \param read_id The read header ID.
     *
     * \param pos The position of the k-mer potentially associated to
     * some variation.
     *
     * \return Return the updated variant information.
     */
    VariantIdentification &add(ReadID_type read_id, size_t pos);

    /**
     * Get the list of read IDs.
     *
     * \return Return the list of the read IDs.
     */
    std::list<ReadID_type> getReads() const;

    /**
     * Get the list of k-mers positions associated to the variant for
     * the given read.
     *
     * \param read_id The read ID.
     *
     * \return Returns the list of positions of the k-mers associated
     * to the variant.
     */
    const std::list<size_t> &getKmersPosInRead(ReadID_type read_id) const;

    /**
     * Get the score of the read associated to the variant.
     *
     * \param read_id The read ID.
     *
     * \return Returns the ratio between the longuest consecutive
     * sequence of k-mer positions in the read and the total number of
     * k-mer positions. If the read id is not found (there is no
     * associated k-mer positions, then returns -1).
     */
    double getReadScore(ReadID_type read_id) const;

    /**
     * Get the bounds of the longuest consecutive sequence from the
     * given list.
     *
     * \param sequence A list of (unique) positions.
     *
     * \return Returns the bounds of the longest consecutive sequence
     * from the given list.
     */
    static std::pair<size_t, size_t> getLonguestSequenceBounds(const std::list<size_t> sequence);

  };

}

#endif
// Local Variables:
// mode:c++
// End:
