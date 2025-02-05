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

#include "monitor.h"

#include "config.h"

#include <cassert>

using namespace std;

BEGIN_KIM_NAMESPACE

string Monitor::abbreviation(MemoryUnit unit) {
  switch (unit) {
  case Byte:      return "B";
  case Kilobyte:  return "KB";
  case Kibibyte:  return "KiB";
  case Megabyte:  return "MB";
  case Mebibyte:  return "MiB";
  case Gigabyte:  return "GB";
  case Gibibyte:  return "GiB";
  case Terabyte:  return "TB";
  case Tebibyte:  return "TiB";
  case Petabyte:  return "PB";
  case Pebibyte:  return "PiB";
  case Exabyte:   return "EB";
  case Exbibyte:  return "EiB";
  default: throw Exception("Memory unit not handled");
  }
}

string Monitor::memoryWithUnit2string(double v, MemoryUnit unit, int precision) {
  string res;
  if ((unit == AutoDecimal) || (unit == AutoBinary)) {
    const bool binary = (unit == AutoBinary);
    const double threshold = (binary ? Kibibyte : Kilobyte);
    unit = Byte;
    while (v > threshold) {
      v /= threshold;
      switch (unit) {
      case Byte: unit = binary ? Kibibyte : Kilobyte; break;
      case Kilobyte: unit = Megabyte; break;
      case Kibibyte: unit = Mebibyte; break;
      case Megabyte: unit = Gigabyte; break;
      case Mebibyte: unit = Gibibyte; break;
      case Gigabyte: unit = Terabyte; break;
      case Gibibyte: unit = Tebibyte; break;
      case Terabyte: unit = Petabyte; break;
      case Tebibyte: unit = Pebibyte; break;
      case Petabyte: unit = Exabyte; break;
      case Pebibyte: unit = Exbibyte; break;
      default: throw Exception("Memory unit not handled");
      }
    }
  } else {
    v /= double(unit);
  }
  char buffer[100] = { 0 };
  sprintf(buffer, "%.*f", precision, v);
  res += buffer;
  res += " ";
  res += abbreviation(unit);
  return res;
}

Monitor::Monitor(int who):
  _rusage_start(), _rusage_stop(),
  _clock_start(), _clock_stop(),
  _paused(0),
  _who(who), _is_running(false) {
  start();
}

Monitor::duration Monitor::getWallClockTime() {
  const bool was_running = isRunning();
  stop();
  const Monitor::duration elapsed_seconds = chrono::duration_cast<std::chrono::seconds>(_clock_stop - _clock_start) - _paused;
  if (was_running) {
    resume();
  }
  return elapsed_seconds;
}

Monitor::duration Monitor::getUserTime() {
  const bool was_running = isRunning();
  stop();
  assert(_rusage_stop.ru_utime.tv_sec >= _rusage_start.ru_utime.tv_sec);
  const Monitor::duration elapsed_seconds
    = 1s * (((_rusage_stop.ru_utime.tv_sec - _rusage_start.ru_utime.tv_sec)
             + ((_rusage_stop.ru_utime.tv_sec >= _rusage_start.ru_utime.tv_sec)
                ? ((_rusage_stop.ru_utime.tv_sec - _rusage_start.ru_utime.tv_sec) >= 500000)
                : -((_rusage_start.ru_utime.tv_sec - _rusage_stop.ru_utime.tv_sec) >= 500000))));
  if (was_running) {
    resume();
  }
  return elapsed_seconds;
}

Monitor::duration Monitor::getSystemTime() {
  const bool was_running = isRunning();
  stop();
  assert(_rusage_stop.ru_stime.tv_sec >= _rusage_start.ru_stime.tv_sec);
  const Monitor::duration elapsed_seconds
    = 1s * ((_rusage_stop.ru_stime.tv_sec - _rusage_start.ru_stime.tv_sec)
            + ((_rusage_stop.ru_stime.tv_sec >= _rusage_start.ru_stime.tv_sec)
               ? ((_rusage_stop.ru_stime.tv_sec - _rusage_start.ru_stime.tv_sec) >= 500000)
               : -((_rusage_start.ru_stime.tv_sec - _rusage_stop.ru_stime.tv_sec) >= 500000)));
  if (was_running) {
    resume();
  }
  return elapsed_seconds;
}

uint64_t Monitor::getMemory() {
  const bool was_running = isRunning();
  stop();
  assert(_rusage_stop.ru_maxrss >= _rusage_start.ru_maxrss);
  const uint64_t mem = _rusage_stop.ru_maxrss - _rusage_start.ru_maxrss;
  if (was_running) {
    resume();
  }
  return mem * Kibibyte;
}


void Monitor::start() {
  _is_running = true;
  _paused = duration::zero();
  _clock_start = clock::now();
  getrusage(_who, &_rusage_start);
}

void Monitor::stop() {
  _clock_stop = clock::now();
  getrusage(_who, &_rusage_stop);
  _is_running = false;
}

void Monitor::resume() {
  if (!_is_running) {
    _paused += chrono::duration_cast<std::chrono::seconds>(clock::now() - _clock_stop);
    _is_running = true;
  }
}

END_KIM_NAMESPACE
