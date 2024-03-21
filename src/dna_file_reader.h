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

#ifndef __DNA_FILE_READER_H__
#define __DNA_FILE_READER_H__

#include <string>
#include <list>
#include <algorithm>

#include <kim_settings.h>
#include <file_reader.h>

namespace kim {

  /**
   * DNA file reader which allows to extract k-mers from sequences.
   */
  class DNAFileReader: public FileReader {

  public:

    /**
     * Available DNA sequence file formats.
     */
    enum Format {

                 /** The Fasta file format defined by Pearson for its eponym software. */
                 FASTA_FORMAT,

                 /** The Fastq file format initially developer by the
                     Wellcome Trust Sanger Institute and as fully
                     describer in the Cock & al article entitled "The
                     Sanger FASTQ file format for sequences with
                     quality scores, and the Solexa/Illumina FASTQ
                     variants", published in Nucleic Acids Res. 2010
                     Apr. 38(6):1767-71 (doi: 10.1093/nar/gkp1137) */
                 FASTQ_FORMAT,

                 /**
                  * The undefined format (an undefined formatted file
                  * is either closed or doesn't allow to extract
                  * k-mers).
                  */
                 UNDEFINED_FORMAT,
    };

    /**
     * The type of some callback function to call each time a new
     * sequence is encountered (see onNewSequence() method).
     */
    typedef void (*onNewSequenceFct)(const DNAFileReader &);

    /**
     * A quite simple structure that stores essential informations to
     * save a file state on some position and allows to restore the
     * file reader to the correct state.
     *
     * \remark Any change in the file may invalidate the informations
     * stored in a FileState variable.
     */
    struct FileState: public FileReader::FileState {

      /**
       * The format of the handled file.
       */
      Format format;

      /**
       * The length of k-mers to extract.
       */
      size_t k;

      /**
       * Flag to distinguish when some new sequence is to be started or
       * if some sequence is currently being processed.
       */
      bool start_symbol_expected;

      /**
       * Number of processed nucleotides in the current sequence in the
       * processed sequence.
       */
      size_t nb_nucleotides;

      /**
       * Number of rightmost consecutive regular (not degenerate)
       * processed nucleotides in the processed sequence.
       */
      size_t nb_consecutive_regular_nucleotides;

      /**
       * The current sequence number in the file (starting from 1).
       */
      size_t current_sequence_id;

      /**
       * The current sequence description.
       */
      std::string current_sequence_description;

      /**
       * The total length of the current sequence (-1 if not known).
       */
      size_t current_sequence_length;

      /**
       * The current k-mer (might be under construction) of the
       * currently processed sequence.
       *
       * For performance reasons this auxiliary k-mer is a circular
       * string starting at _kmer_start_pos.
       */
      std::string kmer_aux;

      /**
       * The starting position of the k-mer in the _kmer string.
       */
      size_t kmer_start_pos;

      /**
       * The current k-mer of the currently processed sequence.
       */
      mutable std::string kmer;

      /**
       * The starting sequence position
       */
      FileReader::FileState sequence_start_file_state;

      /**
       * Default constructor
       */
      FileState();

      /**
       * Reset this DNA file reader state.
       */
      void reset();

    };

  private:

    /**
     * The Current DNA file reader state
     */
    FileState _state;

  protected:

    /**
     * Callback functions to call each time a new sequence is detected
     * (after the header is processed).
     */
    std::list<onNewSequenceFct> _on_new_sequence_cb;


    /**
     * Since this derived class extend the file state properties, the
     * current file state is only accessible through a protected
     * member.
     *
     * \return Return the current file reader state.
     */
    inline virtual FileState &_getState() override {
      return _state;
    }

    /**
     * The hook method called at the end of open().
     *
     * This method check if the first visible character of the file is
     * a sequence starting symbol. On failure, a FileReaderParseError
     * exception is thrown with an explicit message. Be aware that the
     * open() method is called in the base constructor but that in
     * such case, the overriden methods are not known yet thus are not
     * called.
     */
    virtual void _onOpen() override;

    /**
     * The hook method called at the end of close().
     *
     * This method simply set the format to undefined.
     */
    virtual void _onClose() override;

    /**
     * The hook method called at the end of reset().
     *
     * This method overrides the base class one.
     */
    virtual void _onReset() override;

    /**
     * Process the sequence description.
     *
     * The cursor must be located just at the beginning of the new
     * sequence description and the input stream must be valid in
     * order to use this method.
     *
     * This default implementation allows only single line
     * description. If some want to allow multi-lines sequence
     * description (as the original fasta format allows), then some
     * derived class must be defined that overrides this method.
     */
    void _processSequenceDescription();

    /**
     * This method must be called whenever a new sequence is expected.
     *
     * The next visible character is supposed to be one of those given
     * as argument and must be located at the beginning of a new
     * line. If this is not the case, a FileReaderParseError is thrown
     * with an explicit message.
     *
     * \return Returns true if the header has been correctly parsed
     * (i.e., end of file is not reached).
     */
    template<Format>
    bool _parseSequenceDescription();

    /**
     * This method simply properly end the nucleotide sequence reading
     * according to the format (i.e., make the cursor ready to start a
     * new sequence).
     */
    template <Format>
    void _parseEndOfNucleotideSequence();

    /**
     * Parse the file until some new k-mer is available (positioning
     * the cursor at the end of the current k-mer) or some new
     * sequence is being started (positioning the cursor before the
     * sequence starting symbol.
     */
    template<Format>
    void _parse();

    /**
     * The parse function to use according to the file format.
     */
    void (DNAFileReader::*_parse_mth)();

    /**
     * Format specific version that set the reading cursor to the end
     * of the current available sequence (the last position before the
     * new sequence starting character).
     *
     * \param check_consistency When set to false, just goes to the
     * position before the "new sequence character" that starts on a
     * newline without ensuring the validity of the current processed
     * sequence (if any). This performs faster but may lead to bad
     * positioning.
     */
    template <Format>
    void _gotoSequenceEnd(bool check_consistency);

  public:

    /**
     * Creates a DNA sequence file reader.
     *
     * This internally calls the open method.
     *
     * \param k The length of k-mers to extract.
     *
     * \param filename The name of the DNA sequence file to read
     * (calls open() method except if filename is empty).
     *
     * \param warn Emit warnings (or not) while reading the file
     * (default is to emit warnings).
     */
    DNAFileReader(size_t k, const std::string &filename = "", bool warn = true);

    /**
     * Get the detected file format.
     *
     * \return Return the detected file format.
     */
    inline Format getFormat() const {
      return getState().format;
    }

    /**
     * Get the length of extracted k-mers.
     *
     * \return Returns the length of k-mers to extract.
     */
    inline size_t k() const {
      return getState().k;
    }

    /**
     * Get the current k-mer extracted by this file reader object.
     *
     * \return Returns the current available k-mer or an empty string.
     */
    const std::string &getCurrentKmer() const;

    /**
     * Get the k-mer starting at the given relative position in the current sequence.
     *
     * \param p The start position of the expected k-mer in the
     * current sequence (starting from 0). If the position is located
     * before the current position, then restart counting from the
     * beginning of the sequence, which is not optimal at all. If the
     * given positon is the current position, then it returns the
     * current k-mer. If the position is located after the current
     * k-mer position, then the file is read forward to the the
     * desired k-mer. If the given position is greater or equal to the
     * last available k-mer position of the sequence (i.e., the
     * expected k-mer doesn't exist), then an empty string is
     * returned.
     *
     * \param check_consistency When set to false, just goes to the
     * next visible character that should correspond the wanted
     * position, then read the k-mer at this location. If the wanted
     * k-mer doesn't exist in the sequence, the rest of the file
     * reading will be erroneous. When set to true, then it is mostly
     * equivalent to call getNextKmer() until the wanted k-mer is
     * reached (if the wanted position is greater or equal to the
     * sequence length, then the last k-mer is returned. Checking
     * consistency goes slower but is safer.
     *
     * \return Returns the k-mer at the given location or an empty
     * string.
     */
    const std::string &getKmerAt(size_t p, bool check_consistency = true);

    /**
     * Get the next k-mer from this file reader object.
     *
     * This method may update the current sequence description and
     * updates the current k-mer position.
     *
     * \param skip_degenerate When set to true, go to the next k-mer
     * composed of only A, C, G or T/U nucleotides.
     *
     * \return Returns the next available k-mer or an empty string.
     */
    const std::string &getNextKmer(bool skip_degenerate = false);

    /**
     * Check whether the current k-mer contains degenerate symbols.
     *
     * \return Returns true if 1/ the k-mer exists and 2/ doesn't
     * contains degenerate symbols. Otherwise, false is returned.
     */
    bool currentKmerContainsDegenerateSymbols() const {
      return (getState().nb_consecutive_regular_nucleotides < k());
    }

    /**
     * Get the sequence description the current k-mer comes from.
     *
     * \return Returns the current available sequence description or
     * an empty string.
     */
    inline const std::string &getCurrentSequenceDescription() const {
      return getState().current_sequence_description;
    }

    /**
     * Get the current sequence number in the processed file.
     *
     * \return Returns the current sequence number (starting from 1)
     * in the processed file or 0 if no sequence is being processed
     * yet.
     */
    inline size_t getCurrentSequenceID() const {
      return getState().current_sequence_id;
    }

    /**
     * Get the number of read nucleotides in the current sequence.
     *
     * \return Returns the number of processed nucleotides in the
     * current sequence. To obtain the length of the sequence, go to
     * the end of the sequence (gotoSequenceEnd()) then call this
     * method.
     */
    inline size_t getCurrentSequenceProcessedNucleotides() const {
      return getState().nb_nucleotides;
    }

    /**
     * Parse the current sequence until its end in order to compute
     * its length if it is not already known.
     *
     * \return Return the total length of the current sequence.
     */
    size_t computeCurrentSequenceLength();

    /**
     * Get the total length of the current sequence or -1 if not known.
     *
     * \return Return the total length of the current sequence if it
     * is already computed or -1 otherwise.
     */
    inline size_t getCurrentSequenceLength() const {
      return getState().current_sequence_length;
    }

    /**
     * Get the position of the current k-mer in the sequence (starting
     * from 0).
     *
     * \return Returns the position of the current available k-mer or
     * -1.
     */
    inline size_t getCurrentKmerRelativePosition() const {
      return ((getState().nb_nucleotides >= getState().k) ? getState().nb_nucleotides - getState().k : -1);
    }

    /**
     * Set the reading cursor at the end of the current sequence.
     *
     * Calling this method when the cursor is already at the end of
     * the sequence (or on a newly opened file) doesn't change its
     * position.
     *
     * \param check_consistency When set to false, simply goes to the
     * character juste before the next "new sequence character"
     * starting a newline without ensuring the validity of the current
     * processed sequence (if any). This performs faster but may lead
     * to bad positioning.
     */
    void gotoSequenceEnd(bool check_consistency = true);

    /**
     * Set the reading cursor at the beginning of the current
     * sequence.
     *
     * Calling this method when the cursor is already at the beginning
     * of the sequence (or on a newly opened file) doesn't change its
     * position.
     */
    void gotoSequenceStart();

    /**
     * Set the reading cursor to the start of the next available
     * sequence (the first non blank character following the sequence
     * header).
     *
     * Calling this method on a newly file opened file position the
     * cursor to the beginning of the first sequence.
     *
     * \param check_consistency When set to false, just goes to the
     * next "new sequence character" starting a newline without
     * ensuring the validity of the current processed sequence (if
     * any), then read the line as the description and consider the
     * first next non-empty position as the beginning of the
     * sequence. This performs faster but may lead to bad positioning.
     *
     * \return This method returns true if the cursor is correctly
     * positioned to the start of a new sequence and false if no new
     * sequence has been found.
     */
    bool gotoNextSequence(bool check_consistency = true);

    /**
     * Add a callback function to run each time a new sequence is
     * detected (after the description is parsed).
     *
     * \param cb The callback function to call.
     */
    inline void addOnNewSequenceCallback(onNewSequenceFct cb) {
      _on_new_sequence_cb.push_back(cb);
    }

    /**
     * Remove all occurrences of the given callback function (if any)
     * from registered callbacks to run each time a new sequence is
     * detected (after the description is parsed).
     *
     * \param cb The callback function to remove.
     */
    inline void removeOnNewSequenceCallback(onNewSequenceFct cb) {
      std::remove(_on_new_sequence_cb.begin(), _on_new_sequence_cb.end(), cb);
    }

    /**
     * This method is not available for DNA file readers since it
     * would break the current instance state.
     *
     * Thus invoking this method will throw an exception with an
     * explicit message.
     *
     * \return Returns nothing since an exception is always thrown.
     */
    virtual bool setState(const FileReader::FileState &) override final;

    /**
     * Restore this DNA file reader to the given state.
     *
     * If the current file differs from the file defined in the state
     * parameter, then the current file is closed and the file defined
     * in the state parameter is opened.
     *
     * \param s State of the DNA file reader to restore.
     *
     * \return Returns true if the DNA file reader state is correctly
     * restored and false otherwise.
     */
    virtual bool setState(const FileState &s);

    /**
     * Get the current DNA file reader state.
     *
     * \return Returns the current state if the file reader.
     */
    inline virtual const FileState &getState() const override {
      return _state;
    }

  };

  std::ostream &operator<<(std::ostream &os, DNAFileReader::Format fmt);

}

#endif
// Local Variables:
// mode:c++
// End:
