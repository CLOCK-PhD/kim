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

#include <kim_settings.h>

#include <string>
#include <iostream>
#ifdef NDEBUG
#  undef NDEBUG
#endif
#include <cassert>

using namespace std;
using namespace kim;

int main() {

  bool exception_thrown;

  cout << "Testing KimSettings" << endl;

  cout << "Creating default settings" << endl;
  Settings s;
  cerr << "s.alpha = " << s.alpha() << endl;
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(!s.frozen());
  assert(!s.valid());
  assert(s.alpha() == 0.01);
  assert(s.threshold() == 0);
  assert(s.weakMode());
  assert(!s.strictMode());

  cout << "Setting k-mer length to 3" << endl;
  s.setKmerLength(3);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(!s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 3);
  assert(s.alpha() == 0.01);
  assert(s.threshold() == 0);
  assert(s.weakMode());
  assert(!s.strictMode());

  cout << "Setting k-mer prefix length to 1" << endl;
  s.setKmerPrefixLength(1);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(!s.frozen());
  assert(s.getKmerPrefixLength() == 1);
  assert(s.getKmerSuffixLength() == 2);
  assert(s.alpha() == 0.01);
  assert(s.threshold() == 0);
  assert(s.weakMode());
  assert(!s.strictMode());
  assert(s.valid());

  cout << "Trying to set k-mer prefix length to 4 (should lead to 2)" << endl;
  s.setKmerPrefixLength(4);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(!s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 3);
  assert(s.getKmerPrefixLength() == 2);
  assert(s.getKmerSuffixLength() == 1);
  assert(s.alpha() == 0.01);
  assert(s.threshold() == 0);
  assert(s.weakMode());
  assert(!s.strictMode());
  assert(s.valid());

  cout << "Setting k-mer length to 10" << endl;
  s.setKmerLength(10);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(!s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 2);
  assert(s.getKmerSuffixLength() == 8);
  assert(s.alpha() == 0.01);
  assert(s.threshold() == 0);
  assert(s.weakMode());
  assert(!s.strictMode());
  assert(s.valid());

  cout << "Setting k-mer prefix length to 4" << endl;
  s.setKmerPrefixLength(4);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(!s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.01);
  assert(s.threshold() == 0);
  assert(s.weakMode());
  assert(!s.strictMode());
  assert(s.valid());

  cout << "Setting index directory to '/some/path'" << endl;
  s.setIndexDirectory("/some/path");
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory() == "/some/path");
  assert(!s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.01);
  assert(s.threshold() == 0);
  assert(s.weakMode());
  assert(!s.strictMode());
  assert(s.valid());

  for (auto const &p: { "/some/path/that/must/not/exist", "/tmp"}) {
    for (bool must_exist: {false, true}) {
      for (bool must_not_exist: {false, true}) {
        bool expected_exception_thrown = (((p[1] == 's') && must_exist)
                                          || ((p[1] == 't') && must_not_exist));
        exception_thrown = false;
        try {
          cout << "Trying to set index directory to '" << p << "'"
               << " with flag must_exist set to " << must_exist
               << " and flag must_not_exist set to " << must_not_exist << "." << endl
               << (expected_exception_thrown ? "An" : "No")
               << " exception should be thrown..." << endl;
          s.setIndexDirectory(p, must_exist, must_not_exist);
          cout << "Index Directory: '" << s.getIndexDirectory() << "'" << endl;
        } catch (const BadSettingsException &e) {
          cout << "The following exception was thrown: " << e.what() << endl;
          exception_thrown = true;
        }
        if (!exception_thrown) {
          assert(s.getIndexDirectory() == p);
        }
        assert(exception_thrown == expected_exception_thrown);
      }
    }
  }

  cout << "Setting index directory to ''" << endl;
  s.setIndexDirectory("");
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(!s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.01);
  assert(s.threshold() == 0);
  assert(s.weakMode());
  assert(!s.strictMode());
  assert(s.valid());

  cout << "Enable warnings." << endl;
  s.warn(true);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(!s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.01);
  assert(s.threshold() == 0);
  assert(s.weakMode());
  assert(!s.strictMode());
  assert(s.valid());

  cout << "Disable warnings." << endl;
  s.warn(false);
  assert(!s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(!s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.01);
  assert(s.threshold() == 0);
  assert(s.weakMode());
  assert(!s.strictMode());
  assert(s.valid());

  cout << "Enable warnings again." << endl;
  s.warn(true);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(!s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.01);
  assert(s.threshold() == 0);
  assert(s.weakMode());
  assert(!s.strictMode());
  assert(s.valid());

  cout << "Enable consistency checking." << endl;
  s.checkConsistency(true);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(!s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.01);
  assert(s.threshold() == 0);
  assert(s.weakMode());
  assert(!s.strictMode());
  assert(s.valid());

  cout << "Disable consistency checking." << endl;
  s.checkConsistency(false);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(!s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.01);
  assert(s.threshold() == 0);
  assert(s.weakMode());
  assert(!s.strictMode());
  assert(s.valid());

  cout << "Allow overwrite." << endl;
  s.allowOverwrite(true);
  assert(s.warn());
  assert(s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(!s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.01);
  assert(s.threshold() == 0);
  assert(s.weakMode());
  assert(!s.strictMode());
  assert(s.valid());

  cout << "Forbid overwrite." << endl;
  s.allowOverwrite(false);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(!s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.01);
  assert(s.threshold() == 0);
  assert(s.weakMode());
  assert(!s.strictMode());
  assert(s.valid());

  cout << "Set error of type I to 0.05." << endl;
  s.alpha(0.05);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(!s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0);
  assert(s.weakMode());
  assert(!s.strictMode());
  assert(s.valid());

  cout << "Set k-mer rate threshold to 0.15." << endl;
  s.threshold(0.15);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(!s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0.15);
  assert(s.weakMode());
  assert(!s.strictMode());
  assert(s.valid());

  cout << "Set weak mode to false." << endl;
  s.weakMode(false);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(!s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0.15);
  assert(!s.weakMode());
  assert(s.strictMode());
  assert(s.valid());

  cout << "Set strict mode to false." << endl;
  s.strictMode(false);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(!s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0.15);
  assert(s.weakMode());
  assert(!s.strictMode());
  assert(s.valid());

  cout << "Set strict mode to true." << endl;
  s.strictMode(true);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(!s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0.15);
  assert(!s.weakMode());
  assert(s.strictMode());
  assert(s.valid());

  cout << "Freeze settings." << endl;
  s.freeze();
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0.15);
  assert(!s.weakMode());
  assert(s.strictMode());
  assert(s.valid());

  exception_thrown = false;
  try {
    cout << "Trying to disable warnings (should throw a BadSettingsException)" << endl;
    s.warn(false);
  } catch (const BadSettingsException &e) {
    cout << "The following exception was thrown: " << e.what() << endl;
    exception_thrown = true;
  }
  assert(exception_thrown);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0.15);
  assert(!s.weakMode());
  assert(s.strictMode());
  assert(s.valid());

  exception_thrown = false;
  try {
    cout << "Trying to allow file overwite (should throw a BadSettingsException)" << endl;
    s.allowOverwrite(true);
  } catch (const BadSettingsException &e) {
    cout << "The following exception was thrown: " << e.what() << endl;
    exception_thrown = true;
  }
  assert(exception_thrown);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0.15);
  assert(!s.weakMode());
  assert(s.strictMode());
  assert(s.valid());

  exception_thrown = false;
  try {
    cout << "Trying to set k-mer length to 5 (should throw a BadSettingsException)" << endl;
    s.setKmerLength(5);
  } catch (const BadSettingsException &e) {
    cout << "The following exception was thrown: " << e.what() << endl;
    exception_thrown = true;
  }
  assert(exception_thrown);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0.15);
  assert(!s.weakMode());
  assert(s.strictMode());
  assert(s.valid());

  exception_thrown = false;
  try {
    cout << "Trying to set k-mer prefix length to 5 (should throw a BadSettingsException)" << endl;
    s.setKmerPrefixLength(5);
  } catch (const BadSettingsException &e) {
    cout << "The following exception was thrown: " << e.what() << endl;
    exception_thrown = true;
  }
  assert(exception_thrown);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0.15);
  assert(!s.weakMode());
  assert(s.strictMode());
  assert(s.valid());

  exception_thrown = false;
  try {
    cout << "Trying to set Index Directory to '/some/path' (should throw a BadSettingsException)" << endl;
    s.setIndexDirectory("/some/path");
  } catch (const BadSettingsException &e) {
    cout << "The following exception was thrown: " << e.what() << endl;
    exception_thrown = true;
  }
  assert(exception_thrown);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0.15);
  assert(!s.weakMode());
  assert(s.strictMode());
  assert(s.valid());

  exception_thrown = false;
  try {
    cout << "Trying to enable consistency checking (should throw a BadSettingsException)" << endl;
    s.checkConsistency(true);
  } catch (const BadSettingsException &e) {
    cout << "The following exception was thrown: " << e.what() << endl;
    exception_thrown = true;
  }
  assert(exception_thrown);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0.15);
  assert(!s.weakMode());
  assert(s.strictMode());
  assert(s.valid());

  exception_thrown = false;
  try {
    cout << "Trying to disable warnings (should throw a BadSettingsException)" << endl;
    s.warn(false);
  } catch (const BadSettingsException &e) {
    cout << "The following exception was thrown: " << e.what() << endl;
    exception_thrown = true;
  }
  assert(exception_thrown);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0.15);
  assert(!s.weakMode());
  assert(s.strictMode());
  assert(s.valid());

  exception_thrown = false;
  try {
    cout << "Trying to change the error of type I (should throw a BadSettingsException)" << endl;
    s.alpha(0.2);
  } catch (const BadSettingsException &e) {
    cout << "The following exception was thrown: " << e.what() << endl;
    exception_thrown = true;
  }
  assert(exception_thrown);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0.15);
  assert(!s.weakMode());
  assert(s.strictMode());
  assert(s.valid());

  exception_thrown = false;
  try {
    cout << "Trying to change the k-mer rate threshold (should throw a BadSettingsException)" << endl;
    s.threshold(0.2);
  } catch (const BadSettingsException &e) {
    cout << "The following exception was thrown: " << e.what() << endl;
    exception_thrown = true;
  }
  assert(exception_thrown);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0.15);
  assert(!s.weakMode());
  assert(s.strictMode());
  assert(s.valid());

  exception_thrown = false;
  try {
    cout << "Trying to change mode (should throw a BadSettingsException)" << endl;
    s.weakMode(true);
  } catch (const BadSettingsException &e) {
    cout << "The following exception was thrown: " << e.what() << endl;
    exception_thrown = true;
  }
  assert(exception_thrown);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0.15);
  assert(!s.weakMode());
  assert(s.strictMode());
  assert(s.valid());

  cout << "Unfreeze settings." << endl;
  s.unfreeze();
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(!s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(!s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0.15);
  assert(!s.weakMode());
  assert(s.strictMode());
  assert(s.valid());

  cout << "Enable consistency checking." << endl;
  s.checkConsistency(true);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(!s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0.15);
  assert(!s.weakMode());
  assert(s.strictMode());
  assert(s.valid());

  cout << "Freeze settings again." << endl;
  s.freeze();
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0.15);
  assert(!s.weakMode());
  assert(s.strictMode());
  assert(s.valid());

  exception_thrown = false;
  try {
    cout << "Trying to disable warnings (should throw a BadSettingsException)" << endl;
    s.warn(false);
  } catch (const BadSettingsException &e) {
    cout << "The following exception was thrown: " << e.what() << endl;
    exception_thrown = true;
  }
  assert(exception_thrown);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.valid());

  exception_thrown = false;
  try {
    cout << "Trying to set k-mer length to 5 (should throw a BadSettingsException)" << endl;
    s.setKmerLength(5);
  } catch (const BadSettingsException &e) {
    cout << "The following exception was thrown: " << e.what() << endl;
    exception_thrown = true;
  }
  assert(exception_thrown);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0.15);
  assert(!s.weakMode());
  assert(s.strictMode());
  assert(s.valid());

  exception_thrown = false;
  try {
    cout << "Trying to allow file overwite (should throw a BadSettingsException)" << endl;
    s.allowOverwrite(true);
  } catch (const BadSettingsException &e) {
    cout << "The following exception was thrown: " << e.what() << endl;
    exception_thrown = true;
  }
  assert(exception_thrown);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0.15);
  assert(!s.weakMode());
  assert(s.strictMode());
  assert(s.valid());

  exception_thrown = false;
  try {
    cout << "Trying to set k-mer prefix length to 5 (should throw a BadSettingsException)" << endl;
    s.setKmerPrefixLength(5);
  } catch (const BadSettingsException &e) {
    cout << "The following exception was thrown: " << e.what() << endl;
    exception_thrown = true;
  }
  assert(exception_thrown);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0.15);
  assert(!s.weakMode());
  assert(s.strictMode());
  assert(s.valid());

  exception_thrown = false;
  try {
    cout << "Trying to set Index Directory to '/some/path' (should throw a BadSettingsException)" << endl;
    s.setIndexDirectory("/some/path");
  } catch (const BadSettingsException &e) {
    cout << "The following exception was thrown: " << e.what() << endl;
    exception_thrown = true;
  }
  assert(exception_thrown);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0.15);
  assert(!s.weakMode());
  assert(s.strictMode());
  assert(s.valid());

  exception_thrown = false;
  try {
    cout << "Trying to enable consistency checking (should throw a BadSettingsException)" << endl;
    s.checkConsistency(true);
  } catch (const BadSettingsException &e) {
    cout << "The following exception was thrown: " << e.what() << endl;
    exception_thrown = true;
  }
  assert(exception_thrown);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0.15);
  assert(!s.weakMode());
  assert(s.strictMode());
  assert(s.valid());

  exception_thrown = false;
  try {
    cout << "Trying to disable warnings (should throw a BadSettingsException)" << endl;
    s.warn(false);
  } catch (const BadSettingsException &e) {
    cout << "The following exception was thrown: " << e.what() << endl;
    exception_thrown = true;
  }
  assert(exception_thrown);
  assert(s.warn());
  assert(!s.allowOverwrite());
  assert(s.checkConsistency());
  assert(s.getIndexDirectory().empty());
  assert(s.frozen());
  assert(s.getKmerLength() == s.k());
  assert(s.getKmerLength() == 10);
  assert(s.getKmerPrefixLength() == 4);
  assert(s.getKmerSuffixLength() == 6);
  assert(s.alpha() == 0.05);
  assert(s.threshold() == 0.15);
  assert(!s.weakMode());
  assert(s.strictMode());
  assert(s.valid());

  Settings s2(15, 5, "/some/path", true, true, true, 0.2, 0.99999, false, true);
  assert(s2.warn());
  assert(s2.allowOverwrite());
  assert(s2.checkConsistency());
  assert(s2.getIndexDirectory() == "/some/path");
  assert(s2.frozen());
  assert(s2.getKmerLength() == s2.k());
  assert(s2.getKmerLength() == 15);
  assert(s2.getKmerPrefixLength() == 5);
  assert(s2.getKmerSuffixLength() == 10);
  assert(s2.alpha() == 0.2);
  assert(s2.threshold() == 0.99999);
  assert(!s2.weakMode());
  assert(s2.strictMode());
  assert(s2.valid());

  cout << "That's All, Folk!!!" << endl;
  return 0;

}
