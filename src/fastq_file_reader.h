/******************************************************************************
*                                                                             *
*  Copyright © 2023-2024 -- IGH / LIRMM / CNRS / UM                           *
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

#ifndef __FASTQ_FILE_READER_H__
#define __FASTQ_FILE_READER_H__

#include <string>
#include <fstream>

#include <kim_exception.h>
#include <kim_settings.h>
#include <file_reader.h>

namespace kim {

  /**
   * Fastq file reader which allows to extract k-mers from reads.
   */
  class FastqFileReader: public FileReader {

  private:

    /**
     * Flag to distinguish when some new sequence is to be started or
     * if some sequence is currently being processed.
     */
    bool _start_sequence_state;

    /**
     * Number of processed nucleotides in the current sequence in the
     * processed sequence.
     */
    size_t _nb_nucl;

    /**
     * Number of rightmost consecutive valid processed nucleotides in
     * the processed sequence.
     */
    size_t _nb_valid_nucl;

    /**
     * The current sequence header.
     */
    std::string _read_id;

    /**
     * The current k-mer (might be under construction) of the
     * currently processed sequence.
     */
    std::string _kmer;

    /**
     * The hook method called at the end of reset().
     *
     * This method overrides the base class one.
     */
    virtual void _onReset();

    /**
     * Process the sequence header.
     *
     * The cursor must be located just after the '@' symbol starting
     * the new sequence and the input stream must be valid in order to
     * use this method.
     */
    void _processSequenceHeader();

    /**
     * This method must be called whenever a new sequence is expected.
     *
     * The next visible character is supposed to be an '@' located at
     * the beginning of a new line. If this is not the case, a
     * FastqFileReaderParseError is thrown with an explicit message.
     *
     * \return Returns true if the header has been correctly parsed
     * (i.e., end of file is not reached).
     */
    bool _parseSequenceHeader();

    /**
     * Parse the quality separator string and verify its conformity
     * (otherwise a warning is emitted).
     */
    void _parseQualitySeparator();

    /**
     * Simply skip the quality sequence.
     *
     * After this method, either a new sequence is expected or the end
     * of file.
     */
    void _parseQualitySequence();

    /**
     * Parse the file according to its current state.
     *
     * If some new sequence is expected, then process the new sequence
     * header (see _parseSequenceHeader() method). If the quality
     * separator is encountered, the separator is analyzed for
     * conformity (see _parseQualitySeparator() method) then the
     * quality sequence is parsed (see _parseQualitySequence()
     * method). In all other case, the next visible character is
     * expected to be a nucleotide and thus it is processed to update
     * the current k-mer if it is valid and not degenerated.
     */
    void _parse();

  public:

    /**
     * Creates a Fastq file reader with no associated file.
     *
     * \param settings The k-mer identification metric program
     * settings.
     */
    FastqFileReader(const Settings &settings);

    /**
     * Creates a Fastq file reader.
     *
     * This internally calls the open method.
     *
     * \param filename The name of the fastq file to read.
     *
     * \param settings The k-mer identification metric program
     * settings.
     */
    FastqFileReader(const char *filename, const Settings &settings);

    /**
     * Creates a Fastq file reader.
     *
     * This internally calls the open method.
     *
     * \param filename The name of the fastq file to read.
     *
     * \param settings The k-mer identification metric program
     * settings.
     */
    FastqFileReader(const std::string &filename, const Settings &settings);

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
    inline const std::string &getCurrentRead() const {
      return _read_id;
    }

    /**
     * Get the position of the current k-mer in the sequence (starting
     * from 0).
     *
     * \return Returns the position of the current available k-mer or
     * -1.
     */
    inline size_t getCurrentKmerRelativePosition() const {
      return ((_nb_valid_nucl >= _settings.k()) ? _nb_nucl - _settings.k() : -1);
    }

    /**
     * Set the reading cursor to the start of the next available
     * sequence (the first non blank character following the sequence
     * header).
     *
     * Calling this method on a newly file opened file position the
     * cursor to the beginning of the first sequence.
     *
     * \param check_consistency When set to false, just goes to the
     * next '@' character starting a newline without ensuring the
     * validity of the current processed sequence (if any), then read
     * the line as the header and consider the first next non-empty
     * position as the beginning of the sequence. This performs faster
     * but may lead to bad positioning (for example, if some [part of]
     * quality sequence starts with an '@' of some new line).
     *
     * \return This method returns true if the cursor is correctly
     * positioned to the start of a new sequence and false if no new
     * sequence has been found.
     */
    bool gotoNextSequence(bool check_consistency = true);

  };

}

#endif
// Local Variables:
// mode:c++
// End:
