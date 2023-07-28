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

#ifndef __FASTQ_FILE_READER_H__
#define __FASTQ_FILE_READER_H__

//#include <cstdlib>
#include <string>
#include <fstream>

#include <kim_exception.h>
#include <kim_settings.h>

namespace kim {

  class FastqFileReader;

  /**
   * Exception associated to parse error.
   */
  class FastqFileReaderParseError: public Exception {

  public:

    /**
     * Create a parse error exception associated to some fastq file
     * reader.
     *
     * \param reader The fastq file reader associated to this
     * exception.
     */
    FastqFileReaderParseError(FastqFileReader &reader);

  };

  /**
   * Fastq file reader which allows to extract k-mers from reads.
   */
  class FastqFileReader: public std::ifstream {

  private:

    /**
     * k-mer identification metric program settings.
     */
    const Settings &_settings;

    /**
     * The name of the current processed file.
     */
    std::string _filename;

    /**
     * Thus current input stream associated to the current processed
     * file.
     */
    std::ifstream _ifs;

    /**
     * Current line number of the current processed file.
     */
    size_t _line;
    /**
     * Current column number of the current processed file.
     */
    size_t _col;

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
    char _nextVisibleCharacter();

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

    /**
     * The directories to lookup for some filename by default.
     *
     * \see See the findFile() static method.
     */
    static const char *_search_directories[];

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
     * Closes the file stream before destruction.
     */
    ~FastqFileReader();

    /**
     * Get the filename associated to this reader.
     *
     * \return Returns the associated filename.
     */
    inline const std::string &getFilename() const {
      return _filename;
    }

    /**
     * Open the given filename.
     *
     * This method closes the previously associated file stream.
     *
     * \param filename The name of the fastq file to read.
     */
    void open(const std::string &filename);

    /**
     * Close the current input stream if it is opened.
     *
     * This method closes the previously associated file stream.
     */
    void close();

    /**
     * Reset the current reader.
     *
     * This method does nothing if no file is currently being
     * processed. If some file is being processed, this method reset
     * this reader in the same state as it was after opening the file.
     */
    void reset();

    /**
     * Get the current processed line number.
     *
     * \return Return the current line number of the processed file or
     * 0 if no file is opened.
     */
    inline size_t getFileLineNumber() const {
      return _ifs.is_open() ? _line + 1 : 0;
    }

    /**
     * Get the current processed column number.
     *
     * \return Return the current column number of the processed file
     * or 0 if no file is opened.
     */
    inline size_t getFileColumnNumber() const {
      return _ifs.is_open() ? _col + 1 : 0;
    }

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

    /**
     * Cast this reader into boolean.
     *
     * \return This reader is true if its associated input stream is
     * in good() state and is associated to some file.
     */
    inline operator bool() const {
      return _ifs.good() && _ifs.is_open();
    }

    /**
     * Cast this reader into boolean.
     *
     * \return This reader is false if it is... not true.
     */
    inline bool operator !() const {
      return !(bool) *this;
    }

    /**
     * Try to locate the given filename.
     *
     * \param filename The path of the file to locate.
     *
     * \param directories The array of directories to look for the
     * given filename. If not NULL, the array must contain at least
     * one entry and its last entry must be the NULL string (the empty
     * string is not NULL). If not set, the file is searched in the
     * package resource directories.
     *
     * \returns If the given filename exists and is readable, then
     * does nothing except returning its path.  If the filename is not
     * found and is relative (to the current working directory), then
     * it is also searched from the given directories until one file
     * corresponds and is readable.  If no file is found that is
     * readable, an empty string is returned.
     *
     */
    static std::string findFile(const std::string &filename, const char **directories = _search_directories);

  };

}

#endif
// Local Variables:
// mode:c++
// End:
