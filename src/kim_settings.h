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

#ifndef __SETTINGS_H__
#define __SETTINGS_H__

#include <string>

#include <kim_exception.h>

namespace kim {

  /**
   * Bad k-mer identification metric program settings exception.
   */
  class BadSettingsException: public Exception {

  public:

    /**
     * Create an exception dedicated to bad settings.
     *
     * \param msg The initial message of this exception.
     */
    inline BadSettingsException(const std::string &msg = ""): Exception(msg) {}

  };


  /**
   * The k-mer identification metric index settings.
   */
  class Settings {

  private:

    /**
     * The length of the k-mers.
     */
    size_t _k;

    /**
     * The prefix length of the k-mers.
     */
    size_t _p;

    /**
     * The suffix length of the k-mers.
     */
    size_t _s;

    /**
     * The index directory.
     */
    std::string _index_directory;

    /**
     * Emit warning messages or not.
     */
    bool _warn;

    /**
     * Check DNA file consistency (this increase the reading time).
     */
    bool _check_consistency;

    /**
     * Allow file or directory overwriting.
     */
    bool _allow_overwrite;

    /**
     * The type I error (significance) of the variant analysis (the
     * probability to reject a variant that really is in the analyzed
     * data).
     */
    double _alpha;

    /**
     * The threshold used to consider a variant in some read (thus in
     * the whole set of analyzed reads).
     */
    double _threshold;

    /**
     * The settings can't be modified if _frozen is set to true.
     */
    bool _frozen;

  public:

    /**
     * Create a new Settings instance.
     *
     * If k and p have their default value, no verification of
     * settings validity is performed. Thus it is a good practice to
     * ensure validity once settings are corrects (see valid()
     * method).
     *
     * \param k The length of the k-mers.
     *
     * \param p The length of the k-mers prefixes.
     *
     * \param index_directory The index directory.
     *
     * \param warn Activate or deactivate warnings.
     *
     * \param check_consistency Enable or disable DNA file consistency
     * checking while reading.
     *
     * \param allow_overwrite Allow file and directory overwriting.
     *
     * \param alpha The type I error (significance) of the variant
     * analysis (the probability to reject a variant that really is in
     * the analyzed data). The type I error must be in the range [0;
     * 1].
     *
     * \param threshold The k-mer rate threshold to consider a variant
     * being present in some read. The threshold must be in the range
     * [0; 1].
     *
     * \param freeze Do freeze or not the settings.
     */
    Settings(size_t k = 0, size_t p = 0, const std::string &index_directory = "",
             bool warn = true, bool check_consistency = false, bool allow_overwrite = false,
             double alpha = 1., double threshold = 0., bool freeze = false);

    /**
     * Check whether settings are valid or not.
     *
     * \return Returns true if and only if settings are valid.
     */
    inline bool valid() const {
      return ((_p > 0) && (_k > _p) && (_p + _s == _k) && (_alpha >= 0) && (_alpha <= 1) && (_threshold >= 0));
    }

    /**
     * Get the frozen state of current settings.
     *
     * \see See freeze() and unfreeze() methods.
     *
     * \return Returns true if the settings are frozen and false
     * otherwise.
     */
    inline bool frozen() const {
      return _frozen;
    }

    /**
     * Freeze the current settings in order to prevent any further
     * modification.
     *
     * If settings are not valid, then a BadSettingsException exception is thrown.
     *
     * \see See unfreeze(), frozen() and valid() methods.
     */
    void freeze();

    /**
     * Unfreeze the current settings in order to allow modifications.
     *
     * \see See freeze() and frozen() methods.
     */
    inline void unfreeze() {
      _frozen = false;
    }

    /**
     * Get the length of the k-mers (thus the value of k).
     *
     * \remark No verification about settings validity is performed.
     *
     * \see See the k() method.
     *
     * \return Returns the length of the k-mers.
     */
    inline size_t getKmerLength() const {
      return _k;
    }

    /**
     * Get the length of the k-mers (thus the value of k).
     *
     * \remark This is a shortcut for the getKmerLength() method.
     *
     * \return Returns the length of the k-mers.
     */
    inline size_t k() const {
      return getKmerLength();
    }

    /**
     * Set the length of the k-mers (thus the value of k).
     *
     * \remark The given value can't be less than 2 otherwise a
     * BadSettingsException is thrown with an explicit message. If the
     * given value is less or equal to the actual prefix size plus
     * one, then the prefix size is set to k minus one (and a warning
     * is possibly emitted).
     *
     * If the settings are frozen, any attempt to change the length of
     * the k-mers throws a BadSettingsException with an explicit
     * message.
     *
     * \param k The length of the k-mers.
     *
     * \see See warn() methods.
     */
    void setKmerLength(size_t k);

    /**
     * Get the index directory (where files are stored).
     *
     * \remark No verification about settings validity is performed.
     *
     * \return Returns the directory containing the index files.
     */
    inline const std::string &getIndexDirectory() const {
      return _index_directory;
    }

    /**
     * Set the index directory.
     *
     * If the settings are frozen, any attempt to change the index
     * directory throws a BadSettingsException with an explicit
     * message.
     *
     * \param path The directory to use.
     *
     * \param must_exist If true, the index directory must already
     * exists.
     *
     * \param must_not_exist If true, the index directory must not
     * already exists.
     *
     * \see validateDirectory().
     */
    void setIndexDirectory(const std::string &path, bool must_exist = false, bool must_not_exist = false);

    /**
     * Get the prefix length of the k-mers.
     *
     * \remark No verification about settings validity is performed.
     *
     * \return Returns the length of the prefixes of the k-mers.
     */
    inline size_t getKmerPrefixLength() const {
      return _p;
    }

    /**
     * Set the length of the prefixes of the k-mers.
     *
     * \remark The given value can't be zero otherwise a
     * BadSettingsException is thrown with an explicit message. If the
     * given value plus one is greater than k, then the prefix size is
     * set to k minus one (and a warning is possibly emitted).
     *
     * If the settings are frozen, any attempt to change the length of
     * the prefixes of the k-mers throws a BadSettingsException with
     * an explicit message.
     *
     * \param p The length of the prefixes of the k-mers.
     *
     * \see See warn() methods.
     */
    void setKmerPrefixLength(size_t p);

    /**
     * Get the length of the suffixes of the k-mers.
     *
     * This value is simply the difference between the length of the
     * k-mers and the length of the prefixes of the k-mers (the
     * difference is computed when the length or the prefix length of
     * the k-mers is updated to avoid further computations).
     *
     * \remark No verification about settings validity is performed.
     *
     * \return Returns the length of the suffixes of the k-mers.
     */
    inline size_t getKmerSuffixLength() const {
      return _s;
    }

    /**
     * Get the warn state of current settings.
     *
     * \return Returns true if current settings allows warning to be
     * emitted.
     */
    inline bool warn() const {
      return _warn;
    }

    /**
     * Set the warning status for this settings.
     *
     * \param status Warning are enabled only if this parameter is set
     * to true.
     */
    void warn(bool status);

    /**
     * Get the consistency checking state of current settings.
     *
     * \return Returns true if current settings enables consistency
     * checking on DNA file reading.
     */
    inline bool checkConsistency() const {
      return _check_consistency;
    }

    /**
     * Set the consistency checking status for this settings.
     *
     * \param status Consistency checking is enabled only if this
     * parameter is set to true.
     */
    void checkConsistency(bool status);

    /**
     * Get the overwrite allowing state of current settings.
     *
     * \return Returns true if current settings allows file or directory overwriting.
     */
    inline bool allowOverwrite() const {
      return _allow_overwrite;
    }

    /**
     * Set the overwrite allwoing status for this settings.
     *
     * \param status File and directory overwrite is allowed only if
     * this parameter is set to true.
     */
    void allowOverwrite(bool status);

    /**
     * Get the type I error (significance) of the variant analysis
     * (the probability to reject a variant that really is in the
     * analyzed data).
     *
     * \return Returns the type I error.
     */
    inline double alpha() const {
      return _alpha;
    }

    /**
     * Set the type I error (significance) of the variant analysis
     * (the probability to reject a variant that really is in the
     * analyzed data).
     *
     * \param v The type I error (significance) of the variant
     * analysis (the probability to reject a variant that really is in
     * the analyzed data). The type I error must be in the range [0;
     * 1].
     */
    void alpha(double v);

    /**
     * Get the k-mer rate threshold to consider a variant being
     * present in some read.
     *
     * \return Returns the k-mer rate threshold.
     */
    inline double threshold() const {
      return _threshold;
    }

    /**
     * Set the k-mer rate threshold to consider a variant being
     * present in some read.
     *
     * \param v The k-mer rate threshold to consider a variant being
     * present in some read. The threshold must be in the range [0;
     * 1].
     */
    void threshold(double v);

    /**
     * Check whether the given directory meets expected status.
     *
     * \param path The directory to use.
     *
     * \param must_exist If true, the directory must already exists.
     *
     * \param must_not_exist If true, the directory must not already
     * exists.
     *
     * \note Obvisouly, a directory with both must_exist and
     * must_not_exist set to true necessarily leads to a
     * BadSettingsException.
     */
    static void validateDirectory(const std::string &path, bool must_exist = false, bool must_not_exist = false);

  };

}

#endif
// Local Variables:
// mode:c++
// End:
