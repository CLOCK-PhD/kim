/******************************************************************************
*                                                                             *
*  Copyright © 2024      -- IGH / LIRMM / CNRS / UM                           *
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

#ifndef __MONITOR_H__
#define __MONITOR_H__

#include <chrono>
#include <string>
#include <sys/resource.h>
#include <sys/time.h>

#include <kim_exception.h>

namespace kim {

  /**
   * Monitor exception.
   */
  class MonitorException: public Exception {

  public:

    /**
     * Create an exception dedicated to monitoring.
     *
     * \param msg The initial message of this exception.
     */
    inline MonitorException(const std::string &msg = ""): Exception(msg) {}

  };

  /**
   * Helper class to monitor running time and memory usage fo some
   * process.
   */
  class Monitor {

  public:

    /**
     * The duration between two time points.
     */
    typedef std::chrono::duration<long> duration;

    /**
     * The clock to use for monitoring wall clock time.
     */
    typedef std::chrono::steady_clock clock;

    /**
     * The type of time points used by the clock and the duration.
     */
    typedef std::chrono::time_point<clock> time_point;

    /**
     * Handles memory units using both Binary and Decimal SI prefixes.
     */
    enum MemoryUnit: uint64_t {
      AutoDecimal = size_t(-1),    /**< Use the appropriate memory unit using decimal SI prefix */
      AutoBinary = 0,              /**< Use the appropriate memory unit using binary SI prefix */
      Byte = 1,                    /**< Basic unit */
      Kilobyte = Byte * 1000,      /**< \f$1\,\mbox{kilo} = 10^{3}\f$ */
      Kibibyte = Byte << 10,       /**< \f$1\,\mbox{kibi} = 2^{10}\f$ */
      Megabyte = Kilobyte * 1000,  /**< \f$1\,\mbox{mega} = 10^{6} = 10^{3}\,\mbox{kilo}\f$ */
      Mebibyte = Kibibyte << 10,   /**< \f$1\,\mbox{mebi} = 2^{20} = 2^{10}\,\mbox{kibi}\f$ */
      Gigabyte = Megabyte * 1000,  /**< \f$1\,\mbox{giga} = 10^{9} = 10^{3}\,\mbox{mega}\f$ */
      Gibibyte = Mebibyte << 10,   /**< \f$1\,\mbox{gibi} = 2^{20} = 2^{10}\,\mbox{mebi}\f$ */
      Terabyte = Gigabyte * 1000,  /**< \f$1\,\mbox{tera} = 10^{12} = 10^{3}\,\mbox{giga}\f$ */
      Tebibyte = Gibibyte << 10,   /**< \f$1\,\mbox{tebi} = 2^{20} = 2^{10}\,\mbox{gibi}\f$ */
      Petabyte = Terabyte * 1000,  /**< \f$1\,\mbox{peta} = 10^{15} = 10^{3}\,\mbox{tera}\f$ */
      Pebibyte = Tebibyte << 10,   /**< \f$1\,\mbox{pebi} = 2^{20} = 2^{10}\,\mbox{tebi}\f$ */
      Exabyte  = Petabyte * 1000,  /**< \f$1\,\mbox{exa}  = 10^{18} = 10^{3}\,\mbox{peta}\f$ */
      Exbibyte = Pebibyte << 10    /**< \f$1\,\mbox{exbi} = 2^{20} = 2^{10}\,\mbox{pebi}\f$ */
    };

    /**
     * Get the memory unit abbreviation symbol.
     *
     * \param unit The unit for which the abbreviation is wanted.
     *
     * \return Returns the unit abbreviation (following the IEC
     * 60027-2 A.2 and ISO/IEC 80000:13-2008 standards).
     */
    static std::string abbreviation(MemoryUnit unit);

    /**
     * Get a formatted string corresponding to the given value for the
     * given unit, using the given prevision.
     *
     * \param v The value to format (in Bytes).
     *
     * \param unit The unit to use. If the unit is AutoBinary
     * (default) or AutoDecimal, then the unit is choosen in order to
     * have the integer part greater than 0 and less than 1024 for
     * binary SI prefix and 1000 for decimal SI prefix.
     *
     * \param precision The number of decimals in the output string.
     *
     * \return Returns the formatted string corresponding to the given value.
     */
    static std::string memoryWithUnit2string(double v, MemoryUnit unit = AutoBinary, int precision = 3);

  private:

    /**
     * The resource usage on start.
     */
    struct rusage _rusage_start;

    /**
     * The resource usage on stop.
     */
    struct rusage _rusage_stop;

    /**
     * The wall clock time on start.
     */
    time_point _clock_start;

    /**
     * The wall clock time on stop.
     */
    time_point _clock_stop;

    /**
     * The duration of pause(s) since the monitor was (re)started.
     */
    duration _paused;

    /**
     * The process to monitor (see getrusage())
     */
    int _who;

    /**
     * The running state of this monitor.
     */
    bool _is_running;

  public:

    /**
     * Build a monitoring instance.
     *
     * The monitor instance is automatically started.
     *
     * \param who The process to monitor
     *
     * \see getrusage() from the standard library.
     */
    Monitor(int who = RUSAGE_SELF);

    /**
     * Get the wall clock time elapsed since the last time this
     * monitor was started.
     *
     * \return Returns the duration (wall clock) of the monitored
     * process since the last time this monitor was started. All the
     * time the monitor was paused is not taken into account in the
     * result.
     */
    duration getWallClockTime();

    /**
     * Get the user CPU time spent by the monitored process(es).
     *
     * \return Returns the CPU time used by the monitored process(es)
     * in the userland space since the last time this monitor was
     * started.
     */
    duration getUserTime();

    /**
     * Get the system CPU time spent by the monitored process(es).
     *
     * \return Returns the CPU time used by the monitored process(es)
     * in the system space since the last time this monitor was
     * started.
     */
    duration getSystemTime();

    /**
     * Get the memory used (in B) by the monitored process(es).
     *
     * \return Returns the CPU time used by the monitored process(es)
     * in the system space since the last time this monitor was
     * started.
     */
    uint64_t getMemory();

    /**
     * Starts (or restart) the monitoring.
     */
    void start();

    /**
     * Stop (or pause) the monitoring.
     */
    void stop();

    /**
     * Resume the monitoring.
     *
     * The duration between the last stop action and the resume will
     * not be reported by the getWallClockTime() method.
     */
    void resume();

    /**
     * Get this monitor status.
     *
     * \return Returns true if this monitor is running (is not
     * stopped/paused) and false otherwise.
     */
    inline bool isRunning() const {
      return _is_running;
    }

  };

}

#endif
// Local Variables:
// mode:c++
// End:
