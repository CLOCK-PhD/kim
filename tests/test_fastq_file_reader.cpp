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

#include <kim_settings.h>
#include <fastq_file_reader.h>

#include <cstdlib>
#include <iostream>
#include <string>
#include <regex>
#ifdef NDEBUG
#  undef NDEBUG
#endif
#include <cassert>

using namespace std;
using namespace kim;

int main() {

  Settings settings(5, 2);
  assert(settings.valid());
  settings.freeze();
  assert(settings.frozen());
  cout << "Testing FastqFileReader"
       << " with k = " << settings.getKmerLength()
       << " and p = " << settings.getKmerPrefixLength()
       << endl;

  string fname = "test-reads.fastq";

  cout << "Trying to open file '" << fname << "' which is not in the current working directory." << endl;
  FastqFileReader reader(fname, settings);
  cout << "Current opened filename should be empty: '" << reader.getFilename() << "'" << endl;
  assert(reader.getFilename().empty());

  cout << "Looking for file '" << fname << "' in the package directories." << endl;
  const char *dirs[] = {
                      "./",
                      PACKAGE_DATADIR "/",
                      /* The following directories are mostly for development purpose */
                      "resources/",
                      "../resources/",
                      "../",
                      SRCDIR "/../resources/",
                      /* The last array entry must be NULL */
                      NULL
  };
  fname = FastqFileReader::findFile(fname, dirs);
  cout << "File path should not be empty: '" << fname << "'" << endl;
  assert(!fname.empty());

  reader.open(fname);

  cout << "The current filename is '" << reader.getFilename() << "'" << endl;
  assert(reader.getFilename() == fname);

  int cpt = 0;
  cout << "Counting the number of sequences using no checking (thus will be wrong)" << endl;
  while (reader && (cpt < 25)) {
    if (reader.gotoNextSequence(false)) {
      cout << "New sequence: '" << reader.getCurrentRead() << "'" << endl;
      ++cpt;
    } else {
      cout << "No new sequence found." << endl;
    }
  }
  cout << "There is " << cpt << " lines (expecting 21) seen as sequences in '" << reader.getFilename() << "'" << endl;
  assert(cpt == 21);

  cout << "Closing the reader" << endl;
  reader.close();
  assert(!reader);

  vector<size_t> sequence_length;
  sequence_length.reserve(cpt);

  cout << "The current reader associated filename is '" << reader.getFilename() << "' (expecting the empty string)" << endl;
  assert(reader.getFilename().empty());
  cout << "Opening the reader with file '" << fname << "'" << endl;
  reader.open(fname);
  cout << "The current reader associated filename is '" << reader.getFilename() << "' (expecting '" << fname << "')" << endl;
  assert(reader.getFilename() == fname);

  cpt = 0;
  cout << "Counting the number of sequences checking for validity" << endl;
  while (reader) {
    try {
      if (++cpt == 15) {
        cout << "The sequence 15 has a longer quality and should throw an exception." << endl;
      }
      if (reader.gotoNextSequence(true)) {
        cout << "New sequence: '" << reader.getCurrentRead() << "'" << endl;
        // All sequence of the test file have the same pattern:
        regex pattern("Sequence ([0-9]+) of length ([0-9]+) \\(from line ([0-9]+) to line ([0-9]+)\\)(.*)");
        smatch matches;
        assert(regex_match(reader.getCurrentRead(), matches, pattern));
        assert(matches.size() == 6);
        int v = stoi(matches[1]);
        cout << "- Sequence ID is " << v << " (expecting " << cpt << ")" << endl;
        assert(v == cpt);
        v = stoi(matches[2]);
        cout << "- Sequence length is " << v << endl;
        sequence_length.push_back(v);
        v = stoi(matches[3]);
        cout << "- Starting line is " << v << " (expecting " << (reader.getFileLineNumber() - 1) << ")" << endl;
        assert(size_t(v) == (reader.getFileLineNumber() - 1));
      } else {
        cout << "No new sequence found." << endl;
        --cpt;
      }
    } catch (const Exception &e) {
      cout << "An Exception with message '" << e.what() << "' has been thrown" << endl;
      if (cpt == 16) { // The expected exception is thrown when attempting to go to the sequence 16
        --cpt;
      } else {
        return 1;
      }
    }
  }
  cout << "There is " << cpt << " sequences (expecting 18) in '" << reader.getFilename() << "'" << endl;
  assert(cpt == 18);

  cout << "Resetting the reader" << endl;
  reader.reset();
  assert(reader);

  cpt = 0;
  string last_header = "";
  size_t k = settings.k();
  for (string kmer = reader.getNextKmer(); !kmer.empty(); kmer = reader.getNextKmer()) {

    size_t p = reader.getCurrentKmerRelativePosition();
    if (last_header != reader.getCurrentRead()) {
      last_header = reader.getCurrentRead();
      cout << "Starting a new sequence: '" << last_header << "'" << endl;
      ++cpt;
    }
    cout << "Current " << k << "-mer '" << kmer << "'  is at position " << p << endl;
    cout << "(reader.getCurrentKmer() = '" << reader.getCurrentKmer() << "' and sequence length is " << sequence_length[cpt - 1] << ")" << endl;
    assert(kmer == reader.getCurrentKmer());
    assert(p + k <= sequence_length[cpt - 1]);
    switch (cpt) {
    case 1:
    case 2:
    case 3:
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
    case 13:
    case 14:
      break;
    case 4:
      if (((p > 5) && (p < 10)) || ((p > 36) && (p < 50))) {
        cout << "There is a N in this " << k << "-mer. This must not happens." << endl;
        return 1;
      }
      break;
    case 5:
      if (((p > 15) && (p < 25)) || ((p > 35) && (p < 43)) || ((p > 55) && (p < 60))) {
        cout << "There is degeneracy symbols in this " << k << "-mer. This must not happens." << endl;
        return 1;
      }
      break;
    case 15:
      if (p == 25) {
        cout << "This is the last " << k << "-mer of this sequence and the quality sequence is longer than the DNA sequence." << endl
             << "In order to not throw and exception, the method gotoNextSequence(check_consistency: false) is invoked." << endl;
        bool ok = reader.gotoNextSequence(false);
        assert(ok);
      }
      break;
    case 16:
      if (p == 25) {
        cout << "The next sequence contains only N and thus contains no " << k << "-mer." << endl
             << "The next " << k << "-mer should come from the sequence 18." << endl;
        ++cpt;
      }
      break;
    case 18:
      if (p == 0) {
        assert(last_header.substr(9, 2) == "18");
      }
      break;
    default:
      cout << "This situation (cpt == " << cpt << ") should not occur" << endl;
      return 1;
    }
  }

  cout << "That's All, Folk!!!" << endl;
  return 0;

}
