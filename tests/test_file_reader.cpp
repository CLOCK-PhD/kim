/******************************************************************************
*                                                                             *
*  Copyright © 2025      -- IGH / LIRMM / CNRS / UM                           *
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

#include <file_reader.h>

#include <iostream>
#include <string>
#ifdef NDEBUG
#  undef NDEBUG
#endif
#include <cassert>

using namespace std;
using namespace kim;

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

ostream &operator<<(ostream &os, const FileReader::CompressionAlgorithm &a) {
  switch (a) {
  case FileReader::NO_COMPRESSION:
    os << "not compressed";
    break;
  case FileReader::GZIP:
    os << "compressed with gzip";
    break;
  case FileReader::BZ2:
    os << "compressed with bzip2";
    break;
  }
  return os;
}

ostream &operator<<(ostream &os, const FileReader::FileState &s) {
  os << s.filename
     << "[" << s.compression_algorithm << "]"
     << ":" << s.line << "," << s.column
     << "|" << s.pos
     << "(cur: '" << s.current_char << "', next: '" << s.next_char << "')";
  return os;
}

bool operator==(const FileReader::FileState &s1, const FileReader::FileState &s2) {
  return ((s1.filename == s2.filename)
          && (s1.line == s2.line)
          && (s1.column == s2.column)
          && (s1.pos == s2.pos)
          && (s1.compression_algorithm == s2.compression_algorithm)
          && (s1.current_char == s2.current_char)
          && (s1.next_char == s2.next_char));
}


class FileReaderTester {
private:
  FileReader &_reader;
public:
  FileReaderTester(FileReader &reader, const string &fname, FileReader::CompressionAlgorithm expected_compression_algorithm);

  void check(const FileReader::FileState &state) const;
};

FileReaderTester::FileReaderTester(FileReader &reader, const string &fname, FileReader::CompressionAlgorithm expected_compression_algorithm): _reader(reader) {
  cout << "=== " << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__
       << "("
       << "fname = '" << fname << "', "
       << "reader@" << &_reader << ", "
       << "..."
       << ") ==="
       << endl;

  _reader.warn = true; // Turn on warnings

  cout << "  Trying to open file '" << fname << "' which is not in the current working directory." << endl;
  _reader.open(fname);
  cout << "  Current opened filename should be empty: '" << _reader.getFilename() << "'" << endl;
  assert(_reader.getFilename().empty());

  cout << "  Looking for file '" << fname << "' in the package directories." << endl;
  string f = FileReader::findFile(fname, dirs);
  cout << "  File path should not be empty: '" << f << "'" << endl;
  assert(!f.empty());

  _reader.open(f);

  cout << "  The current filename is '" << _reader.getFilename() << "'" << endl;
  assert(_reader.getFilename() == f);

  cout << "  The current file is " << _reader.getState().compression_algorithm << " (expecting " << expected_compression_algorithm << ")" << endl;
  assert(_reader.getState().compression_algorithm == expected_compression_algorithm);

  cout << "  Closing file" << endl;
  _reader.close();
  assert(!_reader);

  cout << "  Opening file again" << endl;
  _reader.open(f);
  assert(_reader);

  cout << endl;
}

void FileReaderTester::check(const FileReader::FileState &expected_state) const {
  cout << "    Checking file state:" << endl
       << "      Reader state is: " << _reader.getState() << endl
       << "      Expected state is: " << expected_state << endl;
  assert(_reader.getState() == expected_state);
}


struct TestFileProps {
  string filename;
  size_t nb_lines;
  size_t nb_bytes;
  FileReader::CompressionAlgorithm compression_algorithm;
};

int main() {

  FileReader reader;

  TestFileProps test_files[] = {
    { "test-variants.vcf",     82, 9140, FileReader::NO_COMPRESSION },
    { "test-variants.vcf.gz",  82, 9140, FileReader::GZIP },
    { "test-reads.fastq",     112, 5622, FileReader::NO_COMPRESSION },
    { "test-reads.fastq.gz",  112, 5622, FileReader::GZIP },
    { "test-reads.fastq.bz2", 112, 5622, FileReader::BZ2}
  };

  for (TestFileProps &p: test_files) {
    FileReaderTester test(reader, p.filename, p.compression_algorithm);
    FileReader::FileState restore_point_state;
    char char_after_restore_point_state = -1;
    FileReader::FileState expected_state;
    size_t nb_lines = 0;
    size_t total = 0;
    expected_state.filename = FileReader::findFile(p.filename, dirs);
    expected_state.line = 0;
    expected_state.column = 0;
    expected_state.pos = 0;
    expected_state.compression_algorithm = p.compression_algorithm;

    test.check(expected_state);

    cout << "\n  Counting lines:" << endl;
    do {
      expected_state.pos += reader.ignore();
      if (reader) {
        ++expected_state.line;
      } else {
        expected_state.column = reader.getState().column;
      }
      expected_state.current_char = reader.get();
      test.check(expected_state);
      if (expected_state.current_char == 10) {
        ++expected_state.line;
      }
    } while (reader);
    cout << "  File has " << expected_state.line << " lines"
         << " (expecting " << p.nb_lines << ")" << endl;
    assert(expected_state.line == p.nb_lines);
    cout << endl;

    nb_lines = expected_state.line;
    expected_state.pos = expected_state.line = expected_state.column = 0;
    cout << "  Reset file" << endl;
    reader.reset();
    test.check(expected_state);

    cout << "\n  Counting bytes:" << endl;
    do {
      string s = reader.getline();
      size_t n = s.length();
      if (reader.good()) {
        expected_state.current_char = 10;
      } else {
        expected_state.current_char = -1;
      }
      cout << "line is '" << s << "' and has size " << n << endl;
      expected_state.column = n;
      expected_state.pos += n + 1;
      test.check(expected_state);
      if (reader.getState().line == 3) {
        restore_point_state = reader.getState();
        char_after_restore_point_state = reader.peek();
        cout << "Saving state: " << restore_point_state
             << " (where next available char is '" << char_after_restore_point_state << "')"
             << endl;
      }
      total += s.length() + 1 /* the end of line */;
      ++expected_state.line;
    } while (reader);
    --total; // The last line has no end of line.
    --expected_state.line; // Same as above
    assert(nb_lines == expected_state.line);
    cout << "File contains " << total << " characters"
         << " (expecting " << p.nb_bytes << ")" << endl;
    assert(total == p.nb_bytes);

    cout << "\n  Restoring file to state " << restore_point_state << endl;
    bool ok = reader.setState(restore_point_state);
    assert(ok);
    test.check(restore_point_state);
    char c = reader.peek();
    cout << "The next char after this point is '" << c << "'"
         << " (expecting '" << char_after_restore_point_state << "')" << endl;
    reader.peek();
    restore_point_state.next_char = char_after_restore_point_state;
    test.check(restore_point_state);
    cout << endl;
  }

  cout << "That's All, Folk!!!" << endl;
  return 0;

}
