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

#ifndef __FILE_READER_H__
#define __FILE_READER_H__

#include <string>
#include <fstream.h>
#include <gzstream.h>
#include <bzstream.h>

#include <kim_exception.h>
#include <kim_settings.h>

namespace kim {

  class FileReader;

  /**
   * Exception associated to parse error.
   */
  class FileReaderParseError: public Exception {

  public:

    /**
     * Create a parse error exception associated to some file reader.
     *
     * \param reader The file reader associated to this exception.
     */
    FileReaderParseError(FileReader &reader);

  };

  /**
   * Generic file reader.
   */
  class FileReader {

  public:

    /**
     * The compression algorithm FileReader instances can read.
     */
    enum CompressionAlgorithm {
      NO_COMPRESSION, /**< File using no compression (plain files). */
      GZIP,           /**< File compressed with gzip. */
      BZ2             /**< File compressed with bzip2. */
    };

    /**
     * A simple structure that stores essential informations to save a
     * file state on some position and allows to restore the file
     * reader to the correct state.
     *
     * \remark Any change in the file may invalidate the informations
     * stored in a FileState variable.
     */
    struct FileState {

      /**
       * The name of the current processed file.
       */
      std::string filename;

      /**
       * The line number corresponding to the pos file offset.
       */
      size_t line;

      /**
       * The column number corresponding to the pos file offset.
       */
      size_t column;

      /**
       * The input stream offset (from the beginning).
       */
      std::ifstream::pos_type pos;

      /**
       * The input stream compression algorithm
       */
      CompressionAlgorithm compression_algorithm;

      /**
       * The ASCII code of the current character (-1 if unknown).
       */
      int current_char;

      /**
       * The ASCII code of the next coming character (-1 if unknown).
       */
      int next_char;

      /**
       * Creates a default FileState.
       *
       * The line, column, current_char, next_char, and pos are set to
       * -1, the compression_algorithm is set to NO_COMPRESSION and
       * filename to the empty string.
       */
      FileState();

    };

    /**
     * Emit warnings (or not) while reading file.
     */
    bool warn;

  private:

    /**
     * The current file state.
     */
    FileState _state;

    /**
     * The current input stream associated to the current processed
     * file if compressed with BZ2.
     */
    compression::ifstream<bz::streambuf> _ifs_bz2;

    /**
     * The current input stream associated to the current processed
     * file if compressed with GZIP.
     */
    compression::ifstream<gz::streambuf> _ifs_gzip;

    /**
     * Thus current input stream associated to the current processed
     * file if not compressed.
     */
    std::ifstream _ifs_none;

    /**
     * Load the next character of file (and updates _state attribute).
     */
    void _loadNextChar();

  protected:

    /**
     * Get the current input stream internal position.
     *
     * \return Returns the position of the (virtual) stream internal
     * cursor (th enumber of bytes from the beginning).
     */
    std::ifstream::pos_type _tellg();

    /**
     * Set the current input stream position.
     *
     * \param p The stream position to set.
     *
     * \return Returns the position of the current input stream.
     */
    std::ifstream::pos_type _seekg(std::ifstream::pos_type p);

    /**
     * Test if this reader is associated to an open file.
     *
     * \return Returns true if exactly one of the input stream is
     * open.
     */
    inline bool _is_open() const {
      return (_ifs_none.is_open() != _ifs_bz2.is_open()) != _ifs_gzip.is_open();
    }

    /**
     * Since some derived class may want to extend the file state
     * properties, the current file state is only accessible through a
     * protected member that can be overriden f necessary.
     *
     * \return Return the current file reader state.
     */
    inline virtual FileState &_getState() {
      return _state;
    }

    /**
     * The hook method called at the end of open().
     *
     * This method may be overriden by derived class if required (be
     * aware that the open() method is called in the base constructor
     * but that in such case, the overriden methods are not known yet
     * thus are not called)
     */
    inline virtual void _onOpen() {}

    /**
     * The hook method called at the end of close().
     *
     * This method may be overriden by derived class if required.
     */
    inline virtual void _onClose() {}

    /**
     * The hook method called at the end of reset().
     *
     * This method may be overriden by derived class if required.
     */
    inline virtual void _onReset() {}

    /**
     * Detect the compression algorithm used for the current input
     * file.
     */
    void _detectCompressionAlgorithm();

    /**
     * Get the next visible character (ASCII code > 32) from the
     * input stream.
     *
     * This method automatically updates the number of lines and columns.
     *
     * If some blank character is found but is different from space,
     * tabulation or newline, a warning may be emitted (according to
     * the settings).
     *
     * \param stop_before If true, the reading cursor is stopped
     * before the next visible characer (and is available through
     * _state.next_char and/or the peek() method). Otherwise, the next
     * visible character is "consumed" (and is available through
     * _state.current_char).
     *
     * \return Returns the next visible character or nul if an error
     * occured or if the end of the stream is reached.
     */
    char _nextVisibleCharacter(bool stop_before = false);

    /**
     * The directories to lookup for some filename by default.
     *
     * \see See the findFile() static method.
     */
    static const char *_search_directories[];

  public:

    /**
     * Creates a file reader with no associated file.
     *
     * \param filename The name of the file to read (calls open()
     * method except if filename is empty).
     *
     * \param warn Emit warnings (or not) while reading the file
     * (default is to emit warnings).
     */
    FileReader(const std::string &filename = "", bool warn = true);

    /**
     * Closes the file stream before destruction.
     */
    virtual ~FileReader();

    /**
     * Get the character at the current reading cursor position.
     *
     * \return Returns the ASCII code of the character or -1 if some
     * error occurs (i.e., end of file).
     */
    int get();

    /**
     * Get the character next to the current reading cursor position.
     *
     * \return Returns the ASCII code of the character or -1 if some
     * error occurs (i.e., end of file).
     */
    int peek();

    /**
     * Get the string from the current reading cursor position and the
     * given delimiter character.
     *
     * \param delim The character used to denote the "end of line".
     *
     * \return Returns the string from the current reading cursor
     * position and the given delimiter character (or the end of
     * file).
     */
    std::string getline(const char delim = '\n');

    /**
     * Skip the content from the current reading cursor position and
     * the given delimiter character.
     *
     * \param delim The character used to denote the "end of line".
     *
     * \return Returns the number of ignored bytes (including the
     * delimiter character).
     */
    size_t ignore(const char delim = '\n');

    /**
     * Get the filename associated to this reader.
     *
     * \return Returns the associated filename.
     */
    inline const std::string &getFilename() const {
      return getState().filename;
    }

    /**
     * Open the given filename.
     *
     * This method closes the previously associated file stream.
     *
     * \param filename The name of the file to read.
     */
    void open(const std::string &filename);

    /**
     * Get this reader status.
     *
     * \return Returns true if the reader has exactly one of its
     * internal stream correctly opened and with a good status.
     */
    bool good() const;

    /**
     * Clear the current input stream.
     */
    void clear();

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
      return _is_open() ? getState().line + 1 : 0;
    }

    /**
     * Get the current processed column number.
     *
     * \return Return the current column number of the processed file
     * or 0 if no file is opened.
     */
    inline size_t getFileColumnNumber() const {
      return _is_open() ? getState().column + 1 : 0;
    }

    /**
     * Restore this file reader to the given state.
     *
     * If the current file differs from the file defined in the state
     * parameter, then the current file is closed and the file defined
     * in the state parameter is opened.
     *
     * \note This may lead to poor performances if input file is
     * compressed.
     *
     * \param s State of the file reader to restore.
     *
     * \return Returns true if the file reader state is correctly
     * restored and false otherwise.
     */
    virtual bool setState(const FileState &s);

    /**
     * Get the current file reader state.
     *
     * \return Returns the current state if the file reader.
     */
    inline virtual const FileState &getState() const {
      return _state;
    }

    /**
     * Cast this reader into boolean.
     *
     * \return This reader is true if its associated input stream is
     * in good() state and is associated to some file.
     */
    inline operator bool() const {
      return good() && _is_open();
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
