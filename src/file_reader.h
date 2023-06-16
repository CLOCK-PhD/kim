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

#ifndef __FILE_READER_H__
#define __FILE_READER_H__

#include <cstdlib>
#include <string>
#include <fstream>

namespace kim {

  /**
   * Fastq file reader which allows to extract k-mers from reads.
   */
  class FileReader: public std::ifstream {

  private:

    std::string filename;
    std::ifstream ifs;
    size_t line, col;
    bool start_sequence_state;
    size_t nb_nucl, nb_valid_nucl;
    size_t k;
    std::string read_id;
    std::string kmer;
    size_t kmer_pos;
    bool warn;

    /**
     * Get the next visible character (ASCII code > 32) from the
     * input stream.
     *
     * This method automatically updates the number of lines and columns.
     *
     * If a some blank character is found but is different from space,
     * tabulation or newline, a warning may be emitted (according to
     * the settings).
     *
     * \return Returns the next visible character or nul if an error
     * occured or if the end of the stream is reached.
     */
    char nextVisibleCharacter();

  public:

    /**
     * Creates a (fastq) file reader.
     *
     * This internally calls the open method.
     *
     * \param filename The name of the fastq file to read.
     *
     * \param k The length of the k-mers to extract.
     *
     * \param warn Activate (default) or deactivate warnings during
     * file processing.
     */
    FileReader(const char *filename, size_t k, bool warn = true);

    /**
     * Closes the file stream before destruction.
     */
    ~FileReader();

    /**
     * Get the filename associated to this reader.
     *
     * \return Returns the associated filename.
     */
    const std::string &getFilename() const;

    /**
     * Open the given filename.
     *
     * This method closes the previously associated file stream.
     *
     * \param filename The name of the fastq file to read.
     *
     * \param k The length of the k-mers to extract.
     *
     * \param warn Activate (default) or deactivate warnings during
     * file processing.
     */
    void open(const char *filename, size_t k, bool warn = true);

    /**
     * Get the current processed line number.
     *
     * \return Return the current line number of the processed file.
     */
    size_t getFileLineNumber() const;

    /**
     * Get the current processed column number.
     *
     * \return Return the current column number of the processed file.
     */
    size_t getFileColumnNumber() const;

    /**
     * Get the next k-mer from this file reader object.
     *
     * This method may update the current read id and updates the
     * current k-mer position.
     *
     * \return Returns the next available k-mer or an empty string.
     */
    const std::string &getNextKmer();

    /**
     * Get the current k-mer extracted by this file reader object.
     *
     * \return Returns the current available k-mer or an empty string.
     */
    const std::string &getCurrentKmer() const;

    /**
     * Get the read ID the current k-mer comes from.
     *
     * \return Returns the current available read ID or an empty
     * string.
     */
    const std::string &getCurrentRead() const;

    /**
     * Get the position of the current k-mer in the sequence (starting
     * from 0).
     *
     * \return Returns the position of the current available k-mer or
     * -1.
     */
    size_t getCurrentKmerRelativePosition() const;

  };

}

#endif
// Local Variables:
// mode:c++
// End:
