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

#include <kim_settings.h>
#include <file_reader.h>

namespace kim {

  /**
   * DNA file reader which allows to extract k-mers from sequences.
   */
  class DNAFileReader: public FileReader {

  protected:

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
     * The current sequence description.
     */
    std::string _current_sequence_description;

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
    virtual void _onReset() override;

    /**
     * Get the available start symbols for a new sequence.
     *
     * Any derived class must override this method.
     *
     * \return Returns the available start symbols for a new sequence.
     */
    virtual const std::string &_sequenceStartSymbols() const = 0;

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
    virtual void _processSequenceDescription();

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
    virtual bool _parseSequenceDescription();

    /**
     * Parse the file until some new k-mer is available (positioning
     * the cursor at the end of the current k-mer) or some new
     * sequence is started (positioning the cursor at the beginning of
     * the nucleotide sequence.
     *
     * Any derived class must override this method.
     */
    virtual void _parse() = 0;

    enum IUPAC {
                IUPAC_UNDEFINED = 0,
                IUPAC_A = 1,
                IUPAC_C = 2,
                IUPAC_G = 4,
                IUPAC_T = 8,
                IUPAC_U = 8,
                IUPAC_GAP = 16,
                IUPAC_R = IUPAC_A | IUPAC_G, // Purine
                IUPAC_Y = IUPAC_C | IUPAC_T, // Pyrimidine
                IUPAC_K = IUPAC_G | IUPAC_T, // Ketones
                IUPAC_M = IUPAC_A | IUPAC_C, // Amino groups
                IUPAC_S = IUPAC_C | IUPAC_G, // Strong
                IUPAC_W = IUPAC_A | IUPAC_T, // Weak
                IUPAC_B = IUPAC_C | IUPAC_G | IUPAC_T, // Not A
                IUPAC_D = IUPAC_A | IUPAC_G | IUPAC_T, // Not C
                IUPAC_H = IUPAC_A | IUPAC_C | IUPAC_T, // Not G
                IUPAC_V = IUPAC_A | IUPAC_C | IUPAC_G, // Not T neither U
                IUPAC_N = IUPAC_A | IUPAC_C | IUPAC_G | IUPAC_T // Any
    };


    /**
     * Get the IUPAC constant associated to the given nucleotide.
     *
     * \param c The nucleotide character (case insensitive).
     *
     * \return Returns the IUPAC constant corresponding to c.
     */
    static IUPAC _toIUPAC(char c);

  public:

    /**
     * Creates a DNA sequence file reader.
     *
     * This internally calls the open method.
     *
     * \param settings The k-mer identification metric program
     * settings.
     *
     * \param filename The name of the DNA sequence file to read
     * (calls open() method except if filename is empty).
     */
    DNAFileReader(const Settings &settings, const std::string &filename = "");

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
     * current sequence (starting from 0). The position must be
     * located at or after the current position. If the given positon
     * is the current position, then it returns the current k-mer. If
     * the position is located before the current k-mer position, then
     * an exception is thrown with an explicit message. If the given
     * position is greater or equal to the last available k-mer
     * position of the sequence (i.e., the expected k-mer doesn't
     * exist), then an empty string is returned.
     *
     * \param check_consistency When set to false, just goes to the
     * next visible character that should correspond the wanted
     * position, then read the k-mer at this location. If the wanted
     * k-mer doesn't exist in the sequence, the rest of the file
     * reading will be erroneous. When set to true, then it is mostly
     * equivalent to call getNextKmer() until the wanted k-mer is
     * reached. This goes slower but can detect inconstant calls by
     * returning an empty k-mer.
     *
     * \return Returns the k-mer located at the given position or an
     * empty string.
     */
    const std::string &getForwardKmer(size_t p, bool check_consistency = true);

    /**
     * Get the next k-mer from this file reader object.
     *
     * This method may update the current sequence description and
     * updates the current k-mer position.
     *
     * \return Returns the next available k-mer or an empty string.
     */
    const std::string &getNextKmer();

    /**
     * Get the sequence description the current k-mer comes from.
     *
     * \return Returns the current available sequence description or
     * an empty string.
     */
    inline const std::string &getCurrentSequenceDescription() const {
      return _current_sequence_description;
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

  };

}

#endif
// Local Variables:
// mode:c++
// End:
