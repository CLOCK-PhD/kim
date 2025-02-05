/******************************************************************************
*                                                                             *
*  Copyright © 2024-2025 -- IGH / LIRMM / CNRS / UM                           *
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

#ifndef __DNA_FILE_INDEX_H__
#define __DNA_FILE_INDEX_H__

#include <iostream>
#include <list>
#include <map>
#include <string>
#include <vector>


#include <dna_file_reader.h>

namespace kim {

  /**
   * Index to retrieve sequences from a collection of DNA Files.
   *
   * Each encountered sequence is supposed to have a **non empty**
   * description and the first *word* of this description is supposed
   * to be the sequence name. Any non empty sequence of non blank
   * characters is considered as a *word*.
   */
  class DNAFileIndex {

  public:

    /**
     * Prints the given index (using an approximative YAML format) to
     * the given output stream.
     *
     * \param os The stream on which to print the index.
     *
     * \param index The index to print.
     *
     * \return Returns the updated stream.
     */
    friend std::ostream &operator<<(std::ostream &os, const DNAFileIndex &index);

    /**
     * Status code of query operations on this index.
     */
    enum Status {
      SUCCESS = 0,              /**< The operation succeed. */
      SEQUENCE_NOT_FOUND,       /**< The sequence is not found. */
      POSITION_NOT_FOUND,       /**< The sequence position is not found. */
      FILE_NOT_FOUND,           /**< The file is not found. */
      FILE_PARSE_ERROR,         /**< The file can't be parsed. */
      SEQUENCE_NAME_DUPLICATED, /**< The sequence name is already indexed. */
    };

  private:

    /**
     * Structure to store the required informations of some location
     * in a DNA file.
     */
    struct _Bookmark {
      /**
       * The input stream offset (from the beginning).
       */
      std::ifstream::pos_type pos;

      /**
       * The line number corresponding to the pos file offset.
       */
      size_t line;

      /**
       * The column number corresponding to the pos file offset.
       */
      size_t column;

      /**
       * The current char (needed to correctly set line and column numbers).
       */
      int current_char;

      /**
       * Builds a bookmark entry.
       *
       * \param state The file state to register.
       */
      _Bookmark(const FileReader::FileState &state);

    };
    friend std::ostream &operator<<(std::ostream &, const _Bookmark &);

    /**
     * Type for representing vector of bookmarks.
     */
    typedef std::vector<_Bookmark> _Bookmarks;

    /**
     * Structure to store the required informations of sequence in a
     * DNA file.
     */
    struct _Sequence {

      /**
       * The current sequence name.
       */
      std::string name;

      /**
       * The current sequence description.
       */
      std::string description;

      /**
       * Number of total nucleotides in the current sequence (-1 if
       * not completely processed).
       */
      size_t nb_nucleotides;

      /**
       * Bookmarks associated to the current sequence (in the
       * ascending order of positions). Each consecutive bookmark are
       * supposed to be separated by exactly _bookmark_distance
       * nucleotides.
       */
      _Bookmarks bookmarks;

      /**
       * Builds a sequence entry.
       *
       * \param description The sequence description.
       *
       * \param name The sequence name. If empty, the name is deduced
       * from the description.
       */
      _Sequence(const std::string &description, const std::string &name = "");

    };

    /**
     * Type for representing vector of sequence informations.
     */
    typedef std::vector<_Sequence> _Sequences;

    /**
     * The association between sequence names (first word of the
     * sequence description) and their location in the sequence
     * collection.
     */
    typedef std::map<std::string, size_t> _SequenceLocations;

    /**
     * Structure to store the required informations of some DNA file.
     */
    struct _File {

      /**
       * The name of the current processed file.
       */
      std::string filename;

      /**
       * The format of the handled file.
       */
      DNAFileReader::Format format;

      /**
       * Number of total nucleotides in the file until last_pos.
       */
      size_t nb_nucleotides;

      /**
       * Sequences in the current file (ordered by ascending ID, thus
       * by their order of appearance in the file).
       */
      _Sequences sequences;

      /**
       * Associations between sequence names and thier positions in
       * the sequences vector.
       */
      _SequenceLocations sequence_locations;

      /**
       * Reading status of the file (true means the file has been
       * entirely processed).
       */
      bool completed;

      /**
       * Builds a file entry.
       *
       * \param state The file state to register.
       *
       * \param state The file state to register.
       */
      _File(const DNAFileReader::FileState &state);

    };

    /**
     * Type for representing vector of files.
     */
    typedef std::vector<_File> _Files;

    /**
     * The association between file names and their location in the
     * file collection.
     */
    typedef std::map<std::string, size_t> _FileLocations;

    /**
     * Type for representing the association between sequence name and
     * the file name they come from.
     */
    typedef std::map<std::string, std::string> _SequenceFileAssociation;

    /**
     * Collection of indexed files
     */
    _Files _files;

    /**
     * Position of files in _files.
     */
    _FileLocations _file_locations;

    /**
     * The association between sequence name and the file the come
     * from.
     */
    _SequenceFileAssociation _sequence_file_association;

    /**
     * Internal reader to parse DNA file.
     */
    mutable DNAFileReader _reader;

    /**
     * File position of current bookmark in _files (-1 if no bookmark
     * loaded).
     */
    size_t _current_file_pos;

    /**
     * Sequence position of current bookmark in _sequences for the
     * current file (-1 if no bookmark loaded).
     */
    size_t _current_sequence_pos;

    /**
     * Bookmark position of current bookmark in _bookmarks for the
     * current sequence in the current file (-1 if no bookmark
     * loaded).
     */
    size_t _current_bookmark_pos;

    /**
     * Get the position of the _File object associated with the given
     * filename.
     *
     * \param filename The filename to retrieve.
     *
     * \return Returns the position of the corresponding _File instance
     * if filename is found or -1 otherwise.
     */
    size_t _getFilePosition(const std::string &filename) const;

    /**
     * Get the position of the _Sequence object associated with the
     * given sequence_name in file instance at position file_pos.
     *
     * \param file_pos The position of the _File instance in the _files
     * attribute.
     *
     * \param sequence_name The name of the sequence to retrieve.
     *
     * \return Returns the position of the corresponding _Sequence
     * instance if sequence_name is found or -1 otherwise. If file_pos
     * is -1 or greater than the _files attribute size, then the
     * sequence is not found.
     */
    size_t _getSequencePosition(const size_t file_pos, const std::string &sequence_name) const;

    /**
     * Get the address of the _Sequence object associated with the
     * sequence at the given position in filename (starting from 0).
     *
     * \param filename The filename having the wanted sequence.
     *
     * \param position The position of the _Sequence instance in the
     * sequences attribute of the _File instance associated to the
     * given filename.
     *
     * \return Returns the address of the corresponding _Sequence
     * instance if it is found or the NULL pointer otherwise.
     */
    const _Sequence *_getSequence(const std::string &filename, size_t position) const;

    /**
     * Get the address of the _Sequence object associated with the
     * sequence at the given position in filename (starting from 0).
     *
     * \param filename The filename having the wanted sequence.
     *
     * \param sequence_name The name of the sequence to retrieve.
     *
     * \return Returns the address of the corresponding _Sequence
     * instance if it is found or the NULL pointer otherwise.
     */
    const _Sequence *_getSequence(const std::string &filename, const std::string &sequence_name) const;

    /**
     * Try to add some new bookmark after the current one.
     *
     * \remark This method expects the internal reader and the
     * _current_file_pos/_current_sequence_pos/_current_bookmark_pos
     * attribute markers to be correctly set and synchronized.
     *
     * \return Returns true if some new bookmark has been correctly
     * added and false otherwise (if the end of sequence -- or file --
     * is reached).
     */
    bool _addNewBookmark();

    /**
     * Add all bookmarks until the end of the current sequence.
     *
     * \remark This method expects the internal reader and the
     * _current_file_pos/_current_sequence_pos/_current_bookmark_pos
     * attribute markers to be correctly set and synchronized.
     */
    void _indexEndOfCurrentSequence();

    /**
     * Synchronizes the internal reader according to the
     * _current_file_pos/_current_sequence_pos/_current_bookmark_pos
     * attributes.
     *
     * \remark This method is declared as const since the internal
     * reader is mutable.
     *
     * \param position If different from -1, don't reset the internal
     * reader if already in the current bookmark section and the
     * position is located after the current one.
     */
    void _syncReaderWithCurrentBookmark(size_t position = (size_t) -1) const;

    /**
     * Restore the internal reader to the last bookmark saved
     * accordingly to the
     * _current_file_pos/_current_sequence_pos/_current_bookmark_pos
     * attributes.
     *
     * If the _current_file_pos attribute is set to -1, then uses the
     * last indexed file.
     *
     * If the _current_sequence_pos attribute is set to -1, then uses
     * the last indexed sequence for the file at _current_file_pos.
     *
     * If the _current_bookmark_pos attribute is set to -1, then uses
     * the last indexed bookmark for the sequence at
     * _current_sequence_pos in file at _current_file_pos.
     *
     * \remark If the
     * _current_file_pos/_current_sequence_pos/_current_bookmark_pos
     * attributes ware already set to values different from -1, this
     * performs almost synchronization of the internal reader (with
     * poorly performances).
     */
    void _loadLastBookmark();

    /**
     * Adds a new sequence entry in the index for the next sequence of
     * the current file handled by the internal reader.
     *
     * \remark This expects the internal reader and the
     * _current_file_pos/_current_sequence_pos/_current_bookmark_pos
     * attributes being synchronized.
     *
     * The new sequence entry has the initial bookmark and
     * the_current_file_pos/_current_sequence_pos/_current_bookmark_pos
     * are updated accordingly.
     *
     * \return Returns the operation status: SUCCESS if the new
     * sequence entry is correctly added, SEQUENCE_NOT_FOUND if no new
     * sequence is available and SEQUENCE_NAME_DUPLICATED if the
     * sequence is already indexed (consider this last case as an
     * error that breaks the index consistency).
     */
    Status _newSequence();

    /**
     * Adds a new file entry in the index for the given filename.
     *
     * If the file exists and was successfully opened, then the
     * _current_file_pos attribute is set to the position of the
     * corresponding _File instance in the _files attribute. Otherwise
     * the _current_file_pos is set to -1 (check the internal reader
     * state to detect whether the file was opened or not). The
     * _current_sequence_pos and _current_sequence_pos attributes are
     * set to -1.
     *
     * \param filename The filename to index.
     *
     * \return Returns the operation status: SUCCESS if the new file
     * entry is correctly added, FILE_NOT_FOUND if no new file is
     * available and FILE_PARSE_ERROR if the file format is not
     * detected (consider the two later cases as errors that breaks
     * the index consistency).
     */
    Status _newFile(const std::string &filename);

    /**
     * Set the internal reader at the given position for the currently
     * bookmarked sequence in the current file.
     *
     * While setting the cursor at the wante position, the index is
     * updated if necessary (missing bookmarks are added for the
     * current sequence).
     *
     * \remark The current file and current sequence informations must
     * be set.
     *
     * This method is intended to be called by the
     * _setInternalReader(const string &, size_t) method.
     *
     * \param position The position in the sequence (starting from 0)
     * to set the internal reader at.
     *
     * \return Returns SUCCESS if the reader is correctly set to the
     * wanted position or POSITION_NOT_FOUND otherwise.
     */
    Status _setInternalReader(size_t position);

    /**
     * Set the internal reader at the given position of the given
     * sequence name in the currently bookmarked file.
     *
     * While setting the cursor at the wanted position of the wanted
     * sequence, the index is updated if necessary (missing bookmarks
     * are added for the missing sequences).
     *
     * \remark The current file information must be set.
     *
     * This method is intended to be called by the
     * _setInternalReader(const string &, const string&, size_t)
     * method and internally calls _setInternalReader(size_t).
     *
     * \param sequence_name The name of the sequence to reach.
     *
     * \param position The position in the sequence (starting from 0)
     * to set the internal reader at.
     *
     * \return Returns SUCCESS if the sequence is found and the
     * internal reader is set to the given position. If the sequence
     * exists but the position is greater than the sequence lenght,
     * the POSITION_NOT_FOUND status is returned. If the sequence is
     * not found, then the SEQUENCE_NOT_FOUND status is returned.
     */
    Status _setInternalReader(const std::string &sequence_name, size_t position);

    /**
     * Set the internal reader at the given position of the given
     * sequence name from the given file.
     *
     * While setting the cursor at the wanted position of the wanted
     * sequence from the given file, the index is updated if necessary
     * (missing bookmarks are added for the missing sequences).
     *
     * This method internally calls _setInternalReader(const string &,
     * size_t)).
     *
     * \param filename The name of the file having the sequence to
     * reach.
     *
     * \param sequence_name The name of the sequence to reach.
     *
     * \param position The position in the sequence (starting from 0)
     * to set the internal reader at.
     *
     * \return Returns SUCCESS if the sequence is found in the given
     * filename and the internal reader is set to the given
     * position. If the sequence exists in the file but the position
     * is greater than the sequence lenght, the POSITION_NOT_FOUND
     * status is returned. If the sequence is not found in the file,
     * then the SEQUENCE_NOT_FOUND status is returned. If the file
     * can't be opened (whatever the reason, including if the format
     * can't be detected), the FILE_NOT_FOUND status is returned. If
     * the file can be opened but contains at least one formatting
     * error, the FILE_PARSE_ERROR status is returned.
     */
    Status _setInternalReader(const std::string &filename, const std::string &sequence_name, size_t position);

    /**
     * Synchronize the given reader with the internal reader.
     *
     * \param reader The reader to synchronize.
     */
    void _setReaderFromInternalReader(DNAFileReader &reader) const;


  public:

    /**
     * Distance (in number of nucleotides) between two bookmarks.
     */
    const size_t bookmark_distance;

    /**
     * Builds an index to allow storing bookmarks over DNA files.
     *
     * This uses an internal DNA file reader of 1-mers.
     *
     * \param bookmark_distance The number of nucleotides between two
     * consecutive bookmarks in the same sequence.
     */
    DNAFileIndex(const size_t bookmark_distance = 1000);

    /**
     * Check file consistency when using the internal reader.
     *
     * By default, the internal reader checks for file consistency
     * since it needs to jump to possibly inexistent k-mer positions
     * (disabling it is at your own risk). This functionnality makes
     * possible to skip a badly formatted sequence (see
     * test_dna_file_index.cpp to see an example).
     *
     * \param check_consistency When true, file consistency is checked
     * while being processed by the internal reader.
     */
    inline void checkConsistency(bool check_consistency) {
      _reader.check_consistency = check_consistency;
    }

    /**
     * Status of file consistency checking when using the internal
     * reader.
     *
     * By default, the internal reader checks for file consistency
     * since it needs to jump to possibly inexistent k-mer positions
     * (disabling it is at your own risk). This functionnality makes
     * possible to skip a badly formatted sequence (see
     * test_dna_file_index.cpp to see an example).
     *
     * \return Returns true if the internal reader checks for file
     * consistency and false otherwise.
     */
    inline bool checkConsistency() const {
      return _reader.check_consistency;
    }

    /**
     * Enable or disable warnings when parsing file with the internal
     * reader.
     *
     * By default, the internal reader emits warnings when parsing
     * files.
     *
     * \param value When true, warnings are emitted by the internal
     * reader while processing files.
     */
    inline void warn(bool value) {
      _reader.warn = value;
    }

    /**
     * Check whether warnings are enabled when parsing file with the
     * internal reader.
     *
     * By default, the internal reader emits warnings when parsing
     * files.
     *
     * \return Returns true if warnings are emitted by the internal
     * reader while processing files and false otherwise.
     */
    inline bool warn() const {
      return _reader.warn;
    }

    /**
     * Synchronize the internal reader (and the index) with the given
     * reader.
     *
     * \param reader The DNA file reader to synchronize with (this may
     * take some time).
     *
     * \return If the internal reader is correctly synchronized, this
     * method returns SUCCESS. If the filename associated to the given
     * reader is not found or if the reader is not opened, the
     * FILE_NOT_FOUND status is returned. If some error occurs while
     * reading the file a FILE_PARSE_ERROR status is returned. Please,
     * consider the two later cases as errors that breaks index
     * consistency.
     */
    Status sync(const DNAFileReader &reader);

    /**
     * Put the reader stream position at the given nucleotide in the
     * given sequence.
     *
     * If the index is up to date (at least to go to the wanted
     * position), then this can goes very faster than navigating with
     * a combination of
     * gotoNextSequence()/gotoSequenceStart()/gotoSequenceEnd()/getKmerAt()
     * calls. If the index is not up to date, it will run a little bit
     * slowly since while reading, the index is dynamically
     * updated. Although, building an index makes no sense if file is
     * read only once and/or sequentially. The main interest is to
     * retrieve efficiently informations when already processed once.
     *
     * \param reader The DNA file reader to update (and to index).
     *
     * \param filename The name of the file containing the wanted
     * sequence.
     *
     * \param sequence_name The name of the wanted sequence.
     *
     * \param position The wanted position in the given sequence
     * (starting from 0).
     *
     * \return If the sequence name exists in the file being processed
     * by the reader and has a total length greater or equal to the
     * wanted position, then the reader reading cursor is updated to
     * the wanted position and this method returns SUCCESS. If the
     * sequence name exists but the given position doesn't, then this
     * returns POSITION_NOT_FOUND. If the sequence name doesn't exist
     * in the file, then this method returns
     * SEQUENCE_NOT_FOUND. Finally, if the filename doesn't exist this
     * method returns FILE_NOT_FOUND.
     */
    Status set(DNAFileReader &reader, const std::string &filename, const std::string &sequence_name, size_t position);

    /**
     * Create bookmarks (if this index is not up to date) for the
     * current entire sequence.
     *
     * \param reader The DNA file reader to update (and to index).
     *
     * \return If the current sequence exists in the file being
     * processed by the reader, then this index is updated (if needed)
     * and this method returns SUCCESS. If this is not possible to
     * process the current sequence, then this method returns
     * SEQUENCE_NOT_FOUND.
     */
    Status indexCurrentSequence(DNAFileReader &reader);

    /**
     * Create bookmarks (if this index is not up to date) for the
     * given file.
     *
     * \param filename The DNA file to index (does almost nothing if
     * the file is already indexed).
     *
     * \return If the given filename is not found, then returns
     * FILE_NOT_FOUND, otherwise, returns SUCCESS.
     */
    Status indexFile(const std::string &filename);

    /**
     * Clears all informations from this index.
     */
    void clear();

    /**
     * Get the number of (potentially partially) indexed files.
     *
     * \return Returns the number of indexed files (even if no
     * bookmark is already computed).
     */
    size_t getNumberOfFiles() const;

    /**
     * Get the set of (potentially partially) indexed file names.
     *
     * \returns Returns the set of indexed file names (even if no
     * bookmark is stored).
     */
    std::list<std::string> getFilenames() const;

    /**
     * Get the number of (potentially partially) indexed sequences in
     * the given filename or globally if the given filename is the
     * empty string.
     *
     * \param filename The filename for which the number if sequences
     * shall be returned. If filename is empty, then the query is done
     * globally over all indexed files.
     *
     * \return Returns the number of sequences in the given indexed
     * file or of all the indexed files if the given filename is
     * empty. If the filename is not empty but is not indexed, returns
     * 0.
     */
    size_t getNumberOfSequences(const std::string &filename = "") const;

    /**
     * Get the set of (potentially partially) indexed sequences in the
     * given filename or globally if the given filename is the empty
     * string.
     *
     * \param filename The filename containing the sequence names that
     * are returned. If filename is empty, then the query is done
     * globally over all indexed files.
     *
     * \return Returns the name of the sequences in the given indexed
     * file or of all the indexed files if the given filename is
     * empty. If the filename is not empty but is not indexed, an
     * empty set is returned.
     */
    std::list<std::string> getSequenceNames(const std::string &filename = "") const;

    /**
     * Get the sequence name of the given sequence ID within the given
     * filename.
     *
     * \param filename The file containing the sequence name to
     * retrieve.
     *
     * \param id The ID of the sequence (starting from 1) to retrieve
     * the name (the first word of the description).
     *
     * \return Returns the name of the sequence having the given ID
     * within the given filename. If the filename is not indexed or if
     * there is no indexed sequence corresponding to this ID, the
     * empty string is returned (notice that any existing sequence
     * should not have an empty name).
     */
    std::string getSequenceName(const std::string &filename, size_t id) const;

    /**
     * Get the sequence description of the given sequence ID within
     * the given filename.
     *
     * \param filename The file containing the sequence description to
     * retrieve.
     *
     * \param id The ID of the sequence (starting from 1) to retrieve
     * the description.
     *
     * \return Returns the description of the sequence having the
     * given ID within the given filename. If the filename is not
     * indexed or if there is no indexed sequence corresponding to
     * this ID, the empty string is returned (notice that the sequence
     * name should never be the empty string).
     */
    std::string getSequenceDescription(const std::string &filename, size_t id) const;

    /**
     * Get the description of the given sequence within the given
     * filename.
     *
     * \param filename The file containing the description to
     * retrieve.
     *
     * \param sequence_name The name of the sequence to retrieve the
     * description.
     *
     * \return Returns the description of the sequence having the
     * given name within the given filename. If the filename is not
     * indexed or if there is no indexed sequence corresponding to
     * this name, the empty string is returned (notice that the
     * sequence name should never be the empty string).
     */
    std::string getSequenceDescription(const std::string &filename, const std::string &sequence_name) const;

    /**
     * Get the description of the given sequence.
     *
     * This internally call the version with the explicit filename,
     * but needs to first retrieve the filename, thus it takes an
     * extra step of computation.
     *
     * \param sequence_name The name of the sequence to retrieve the
     * description.
     *
     * \return Returns the description of the sequence having the
     * given name. If there is no indexed sequence corresponding to
     * the given name, the empty string is returned (notice that the
     * sequence name should never be the empty string).
     */
    std::string getSequenceDescription(const std::string &sequence_name) const;

    /**
     * Get the ID of the given sequence within the given filename.
     *
     * \param filename The file containing the ID to retrieve.
     *
     * \param sequence_name The name of the sequence to retrieve the
     * ID.
     *
     * \return Returns the ID (starting from 1) of the sequence having
     * the given name within the given filename. If the filename is
     * not indexed or if there is no indexed sequence corresponding to
     * this name, this returns 0.
     */
    size_t getSequenceID(const std::string &filename, const std::string &sequence_name) const;

    /**
     * Get the ID of the given sequence.
     *
     * This internally call the version with the explicit filename,
     * but needs to first retrieve the filename, thus it takes an
     * extra step of computation.
     *
     * \param sequence_name The name of the sequence to retrieve the
     * ID.
     *
     * \return Returns the ID (starting from 1) of the sequence having
     * the given name (the ID is file specific). If there is no
     * indexed sequence corresponding to this name, this returns 0.
     */
    size_t getSequenceID(const std::string &sequence_name) const;

    /**
     * Get the file that contains the given sequence name.
     *
     * \param sequence_name The sequence sequence_name.
     *
     * \return Returns the filename that contains the given sequence
     * name. If there is no indexed sequence corresponding to this
     * name, the empty string is returned.
     */
    std::string getSequenceFile(const std::string &sequence_name) const;

    /**
     * Get the total number of nucleotides processed by this index.
     *
     * \return Returns the number of nucleotides processed by this
     * index.
     */
    size_t getNumberOfNucleotides() const;

    /**
     * Get the total number of nucleotides processed by this index for
     * the given filename.
     *
     * \param filename The file for which the number of nucleotides
     * must be reported.
     *
     * \return Returns the number of nucleotides processed by this
     * index for the given filename.
     */
    size_t getNumberOfNucleotides(const std::string &filename) const;

    /**
     * Get the total number of nucleotides processed by this index for
     * the given sequence.
     *
     * \param filename The file for which the number of nucleotides
     * must be reported.
     *
     * \param sequence_name The sequence name for which the number of
     * nucleotides must be reported.
     *
     * \return Returns the number of nucleotides processed by this
     * index for the given sequence (0 if file or sequence is not found).
     */
    size_t getNumberOfNucleotides(const std::string &filename, const std::string &sequence_name) const;

    /**
     * Get the format of the given file.
     *
     * \param filename The file for which the format is retrieved.
     *
     * \return Returns the format of the given file. If the filename
     * is not indexed, then the DNAFileReader::UNDEFINED_FORMAT is
     * returned.
     */
    DNAFileReader::Format getFileFormat(const std::string &filename) const;

    /**
     * Get the name from the sequence description.
     *
     * The name is actually the first word of the description, knowing
     * that any non empty sequence of non blank characters is
     * considered as a *word*.
     *
     * \param description The sequence description.
     *
     * \return Returns the first word of the description. If the
     * description is empty, then an empty name is returned, but this
     * may lead to index inconsistency.
     */
    static std::string getNameFromDescription(const std::string &description);

    /**
     * Return the string version of the given status.
     *
     * \param status The status code.
     *
     * \return Returns the string representation of the given status.
     */
    static std::string statusString(Status status);

  };

}

#endif
// Local Variables:
// mode:c++
// End:
