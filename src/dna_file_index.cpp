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

#include "dna_file_index.h"

#include "config.h"

#include <cassert>
#include <iostream>

using namespace std;

BEGIN_KIM_NAMESPACE

ostream &operator<<(ostream &os, const FileReader::FileState &s) {
  os << "'" << s.filename << ":" << s.line << ":" << s.column;
  return os;
}

ostream &operator<<(ostream &os, const DNAFileReader::FileState &s) {
  os << "DNAFile[" << s.format << "]:" << (FileReader::FileState) s << "\n"
     << "- k: " << s.k
     << ", start_symbol_expected: " << s.start_symbol_expected
     << ", nb_nucleotides: " << s.nb_nucleotides
     << ", current_sequence_id: " << s.current_sequence_id
     << ", current_sequence_length: " << s.current_sequence_length
     << ", kmer_aux: " << s.kmer_aux
     << ", kmer_start_pos: " << s.kmer_start_pos
     << ", kmer: " << s.kmer
     << ", sequence_start_file_state: " << s.sequence_start_file_state
     << "\n";
  return os;
}

ostream &operator<<(ostream &os, const DNAFileIndex::_Bookmark &b) {
  os << "Bookmark:" << b.line << ":" << b.column << ":" << b.pos;
  return os;
}

ostream &operator<<(ostream &os, const DNAFileIndex &index) {
  os << "Index:";
  if (index._files.empty()) {
    os << " <Empty>" << endl;
  } else {
    os << endl;
    for (size_t f = 0; f < index._files.size(); ++f) {
      const DNAFileIndex::_File &file = index._files[f];
      os << "- Filename: '" << file.filename << "'" << endl
         << "  Format: " << file.format << endl
         << "  Start position: " << file.start_pos << endl
         << "  Last position: " << file.last_pos << endl
         << "  Number of nucleotides: " << file.nb_nucleotides << endl
         << "  Completed: " << (file.completed ? "true" : "false") << endl
         << "  Sequences:" << endl;
      const DNAFileIndex::_Sequences &sequences = file.sequences;
      for (size_t s = 0; s < sequences.size(); ++s) {
        const DNAFileIndex::_Sequence &sequence = sequences[s];
        os << "  - Name: '" << sequence.name << "':" << endl
           << "    Description: '" << sequence.description << "'" << endl
           << "    Length: " << sequence.nb_nucleotides << endl
           << "    Bookmarks: " << endl;
        const DNAFileIndex::_Bookmarks &bookmarks = sequence.bookmarks;
        size_t nb = 0;
        for (size_t b = 0; b < bookmarks.size(); ++b) {
          const DNAFileIndex::_Bookmark &bookmark = bookmarks[b];
          os << "    - " << nb << ": (" << bookmark.line << ", " << bookmark.column << ")" << endl;
          nb += index.bookmark_distance;
        }
      }
    }
  }
  return os;
}

#ifndef NDEBUG
# define ASSERT_VALID_CURRENT_FILE(check_reader_sync) do {              \
    assert(!_files.empty());                                            \
    assert(_current_file_pos != (size_t) -1);                           \
    assert(_current_file_pos < _files.size());                          \
    if (check_reader_sync) {                                            \
      assert(_reader);                                                  \
      const _File &file = _files[_current_file_pos];                    \
      assert(file.filename == _reader.getFilename());                   \
    }                                                                   \
  } while (0)
# define ASSERT_VALID_CURRENT_SEQUENCE(check_reader_sync) do {          \
    ASSERT_VALID_CURRENT_FILE(check_reader_sync);                       \
    const _File &file = _files[_current_file_pos];                      \
    const _Sequences &sequences = file.sequences;                       \
    assert(!sequences.empty());                                         \
    assert(_current_sequence_pos != (size_t) -1);                       \
    assert(_current_sequence_pos < sequences.size());                   \
    if (check_reader_sync) {                                            \
      const _Sequence &sequence = file.sequences[_current_sequence_pos]; \
      assert(sequence.description == _reader.getCurrentSequenceDescription()); \
    }                                                                   \
  } while (0)
# define ASSERT_VALID_CURRENT_BOOKMARK(check_reader_sync) do {          \
    ASSERT_VALID_CURRENT_SEQUENCE(check_reader_sync);                   \
    const _Bookmarks &bookmarks = _files[_current_file_pos].sequences[_current_sequence_pos].bookmarks; \
    assert(!bookmarks.empty());                                         \
    assert(_current_bookmark_pos != (size_t) -1);                       \
    assert(_current_bookmark_pos < bookmarks.size());                   \
    if (check_reader_sync) {                                            \
      size_t nb = (bookmarks.size() + 1) * bookmark_distance;           \
      assert(_reader.getCurrentSequenceProcessedNucleotides() < nb);    \
    }                                                                   \
  } while (0)
#else
# define ASSERT_VALID_CURRENT_FILE(check_reader_sync)   \
  (void) 0
# define ASSERT_VALID_CURRENT_SEQUENCE(check_reader_sync)       \
  (void) 0
# define ASSERT_VALID_CURRENT_BOOKMARK(check_reader_sync)       \
  (void) 0
#endif

DNAFileIndex::_Bookmark::_Bookmark(const FileReader::FileState &state):
  pos(state.pos), line(state.line), column(state.column) {
}

DNAFileIndex::_Sequence::_Sequence(const string &description, const string &name):
  name(name.empty() ? getNameFromDescription(description) : name),
  description(description), nb_nucleotides(-1), bookmarks() {
}

DNAFileIndex::_File::_File(const DNAFileReader::FileState &state):
  filename(state.filename), format(state.format),
  start_pos(state.sequence_start_file_state), last_pos(state.sequence_start_file_state),
  nb_nucleotides(0),
  sequences(), sequence_locations(),
  completed(false) {
}

size_t DNAFileIndex::_getFilePosition(const string &filename) const {
  size_t p = -1;
  // Lookup for the given filename
  const _FileLocations::const_iterator it = _file_locations.find(filename);
  if (it != _file_locations.end()) {
    // The filename exists and is at position it->second in _files
    assert(it->second < _files.size());
    assert(it->first == filename);
    assert(_files[it->second].filename == filename);
    p = it->second;
  }
  return p;
}

size_t DNAFileIndex::_getSequencePosition(const size_t file_pos, const string &sequence_name) const {
  size_t p = -1;
  if (file_pos != (size_t) -1) {
    assert(file_pos < _files.size());
    const _SequenceLocations &seq_loc = _files[file_pos].sequence_locations;
    // Lookup for the given sequence name
    const _SequenceLocations::const_iterator it = seq_loc.find(sequence_name);
    if (it != seq_loc.end()) {
      // The sequence name exists and is at position it->second in seq_loc
      assert(it->second < _files[file_pos].sequences.size());
      assert(it->first == sequence_name);
      assert(_files[file_pos].sequences[it->second].name == sequence_name);
      p = it->second;
    }
  }
  return p;
}

const DNAFileIndex::_Sequence *DNAFileIndex::_getSequence(const string &filename, size_t position) const {
  const _Sequence *ptr = NULL;
  size_t p = _getFilePosition(filename);
  if (p != (size_t) -1) {
    assert(p < _files.size());
    if (position < _files[p].sequences.size()) {
      ptr = &(_files[p].sequences[position]);
    }
  }
  return ptr;
}

const DNAFileIndex::_Sequence *DNAFileIndex::_getSequence(const string &filename, const string &sequence_name) const {
  const _Sequence *ptr = NULL;
  size_t p = _getFilePosition(filename);
  size_t position = _getSequencePosition(p, sequence_name);
  if (position != (size_t) -1) {
    assert(p < _files.size());
    assert(position < _files[p].sequences.size());
    ptr = &(_files[p].sequences[position]);
  }
  return ptr;
}

bool DNAFileIndex::_addNewBookmark() {
  ASSERT_VALID_CURRENT_BOOKMARK(true);
  assert(!_files[_current_file_pos].completed);
  _File &file = _files[_current_file_pos];
  _Sequence &sequence = file.sequences[_current_sequence_pos];
  if (sequence.nb_nucleotides == (size_t) -1) {
    // The sequence lenght is not known, thus there might be some
    // bookmark(s) to add.
    _Bookmarks &bookmarks = sequence.bookmarks;
    size_t next_pos = bookmarks.size() * bookmark_distance;
    assert(next_pos >= bookmark_distance);
    size_t cur_pos = next_pos - bookmark_distance;
    const DNAFileReader::FileState &state = _reader.getState();
    assert(_reader.getCurrentSequenceProcessedNucleotides() >= cur_pos);
    assert(_reader.getCurrentSequenceProcessedNucleotides() < next_pos);
    // Since the internal reader reads 1-mers, the bookmark for position
    // next_pos corresponds to the k-mer at position next_pos - 1.
    if (_reader.getKmerAt(next_pos - 1).empty()) {
      assert(next_pos >= bookmark_distance);
      assert(_reader.getCurrentSequenceLength() != (size_t) -1);
      assert(_reader.getCurrentSequenceLength() >= cur_pos);
      assert(_reader.getCurrentSequenceLength() < next_pos);
      sequence.nb_nucleotides = _reader.getCurrentSequenceLength();
      file.nb_nucleotides += _reader.getCurrentSequenceLength() - cur_pos;
      if (!_reader) {
        file.completed = true;
        _reader.reset();
        _loadLastBookmark();
        _reader.gotoSequenceEnd();
      }
    } else {
      bookmarks.emplace_back(state);
      file.nb_nucleotides += bookmark_distance;
      ++_current_bookmark_pos;
      assert(_current_bookmark_pos + 1 == bookmarks.size());
    }
    file.last_pos = state;
#ifndef NDEBUG
  } else {
    // The sequence lenght is known, thus there shouldn't have any new bookmark to add.
    _Bookmarks &bookmarks = sequence.bookmarks;
    size_t next_pos = bookmarks.size() * bookmark_distance;
    assert(next_pos >= bookmark_distance);
    size_t cur_pos = next_pos - bookmark_distance;
    assert(cur_pos <= sequence.nb_nucleotides);
    assert(next_pos > sequence.nb_nucleotides);
#endif
  }
  return (sequence.nb_nucleotides == (size_t) -1);
}

void DNAFileIndex::_indexEndOfCurrentSequence() {
  while (_addNewBookmark());
}


void DNAFileIndex::_syncReaderWithCurrentBookmark() const {
  ASSERT_VALID_CURRENT_BOOKMARK(false);

  const _File &file = _files[_current_file_pos];
  const _Sequence &sequence = file.sequences[_current_sequence_pos];
  const _Bookmarks &bookmarks = sequence.bookmarks;
  const _Bookmark &first_bookmark = bookmarks.front();
  const _Bookmark &bookmark = bookmarks[_current_bookmark_pos];

  if (_reader.getFilename() != file.filename) {
    _reader.open(file.filename);
  }
  assert(_reader);

  // Retrieve current internal reader state.
  DNAFileReader::FileState state = _reader.getState();

  // Since a bookmark exists, this is not the beginning of the file
  state.start_symbol_expected = false;

  // Informations derived from current sequence
  state.nb_nucleotides = _current_bookmark_pos * bookmark_distance;
  state.nb_consecutive_regular_nucleotides = 0;
  state.current_sequence_length = sequence.nb_nucleotides;
  state.current_sequence_id = _current_sequence_pos + 1;
  state.current_sequence_description = sequence.description;

  // Informations derived from first bookmark of the current sequence
  state.sequence_start_file_state.pos = first_bookmark.pos;
  state.sequence_start_file_state.line = first_bookmark.line;
  state.sequence_start_file_state.column = first_bookmark.column;

  // Informations derived from current bookmark of the current sequence
  state.pos = bookmark.pos;
  state.line = bookmark.line;
  state.column = bookmark.column;

  // Information that can't be retrieved.
  state.kmer = state.kmer_aux = string(state.k, '?');

  // Now set the reader to the wanted bookmark.
  _reader.setState(state);
    assert(_reader.getState().pos != -1); // TEMP

}

void DNAFileIndex::_loadLastBookmark() {

  assert(!_files.empty());
  if (_current_file_pos == (size_t) -1) {
    assert(_current_sequence_pos == (size_t) -1);
    assert(_current_bookmark_pos == (size_t) -1);
    _current_file_pos = _files.size() - 1;
  }

  const _Sequences &sequences = _files[_current_file_pos].sequences;

  assert(!sequences.empty());
  if (_current_sequence_pos == (size_t) -1) {
    assert(_current_bookmark_pos == (size_t) -1);
    _current_sequence_pos = sequences.size() - 1;
  }

  const _Bookmarks &bookmarks = sequences[_current_sequence_pos].bookmarks;

  assert(!bookmarks.empty());
  if (_current_bookmark_pos == (size_t) -1) {
    _current_bookmark_pos = bookmarks.size() - 1;
  }

  _syncReaderWithCurrentBookmark();

}

DNAFileIndex::Status DNAFileIndex::_newSequence() {
  DNAFileIndex::Status status = _reader ? SUCCESS : SEQUENCE_NOT_FOUND;

  if (status == SUCCESS) {

    ASSERT_VALID_CURRENT_FILE(true);

    _File &file = _files[_current_file_pos];
    assert(!file.completed);

    if (_reader.gotoNextSequence()) {
      assert(_reader.getState().pos != -1); // TEMP

      const DNAFileReader::FileState &state = _reader.getState();

      const string &desc = _reader.getCurrentSequenceDescription();
      const string &name = getNameFromDescription(desc);

      if (_sequence_file_association.find(name) == _sequence_file_association.end()) {

        _sequence_file_association[name] = file.filename;

        assert(file.sequences.size() == file.sequence_locations.size());
        assert(file.sequence_locations.find(name) == file.sequence_locations.end());

        file.sequences.emplace_back(desc, name);
        // The inequality is for the case there are some sequence name
        // that are already defined in the index.
        assert(file.sequences.size() <= _reader.getCurrentSequenceID());
        _current_sequence_pos = file.sequences.size() - 1;
        file.sequence_locations[name] = _current_sequence_pos;

        assert(_reader.getCurrentSequenceProcessedNucleotides() == 0);
        file.sequences.back().bookmarks.emplace_back(state);
        _current_bookmark_pos = 0;

      } else {
        status = SEQUENCE_NAME_DUPLICATED;
        _current_sequence_pos = _current_bookmark_pos = (size_t) -1;
      }
      file.last_pos = state;
    } else {
      file.completed = true;
      status = SEQUENCE_NOT_FOUND;
    }
  }

  return status;

}

DNAFileIndex::Status DNAFileIndex::_newFile(const string &filename) {
  DNAFileIndex::Status status = SUCCESS;
  assert(_file_locations.find(filename) == _file_locations.end());
  assert(_files.size() == _file_locations.size());
  _current_file_pos = (size_t) -1;
  _current_sequence_pos = (size_t) -1;
  _current_bookmark_pos = (size_t) -1;
  _reader.open(filename);
  if (_reader) {
    const DNAFileReader::FileState &state = _reader.getState();
    if (state.format != DNAFileReader::UNDEFINED_FORMAT) {
      _File f(state);
      assert(f.sequences.empty());
      assert(f.sequence_locations.empty());
      _files.push_back(f);
      assert(_files.size() > 0);
      assert(_file_locations.find(filename) == _file_locations.end());
      _file_locations[f.filename] = _files.size() - 1;
      assert(_file_locations.find(filename) != _file_locations.end());
      assert(_files.size() == _file_locations.size());
      _current_file_pos = _files.size() - 1;
    } else {
      status = FILE_PARSE_ERROR;
    }
  } else {
    status = FILE_NOT_FOUND;
  }
  return status;
}

DNAFileIndex::Status DNAFileIndex::_setInternalReader(size_t position) {
  ASSERT_VALID_CURRENT_SEQUENCE(false);
  _Sequence &sequence = _files[_current_file_pos].sequences[_current_sequence_pos];
  const _Bookmarks &bookmarks = sequence.bookmarks;
  assert(!bookmarks.empty());

  size_t expected_bookmark_pos = (position
                                  ? ((position - 1) / bookmark_distance)
                                  : 0);

  if (sequence.nb_nucleotides != (size_t) -1) {
    // The sequence is already fully indexed
    if (position > sequence.nb_nucleotides) {
      // The position doesn't exist, Setting the internal reader to
      // the end of the current sequence.
      _current_bookmark_pos = bookmarks.size() - 1;
      _syncReaderWithCurrentBookmark();
    } else {
      // The position exists, thus going to the closest bookmark of
      // the current sequence.
      assert(expected_bookmark_pos <= bookmarks.size());
      _current_bookmark_pos = expected_bookmark_pos;
      _syncReaderWithCurrentBookmark();
    }
  } else {
    // The sequence is not fully indexed.
    //
    // Going to the last bookmark for this sequence
    _current_bookmark_pos = bookmarks.size() - 1;
    _syncReaderWithCurrentBookmark();
    // and update missing bookmarks if necessary.
    while ((_current_bookmark_pos < expected_bookmark_pos) && _addNewBookmark());
  }

  Status status = SUCCESS;
  ASSERT_VALID_CURRENT_SEQUENCE(true);

  if (_current_bookmark_pos < expected_bookmark_pos) {
    // The whole sequence has been indexed, but the position is too
    // far away the last nucleotide.
    assert(sequence.nb_nucleotides != (size_t) -1);
    assert(_reader.getCurrentSequenceLength() != (size_t) -1);
    assert(sequence.nb_nucleotides == _reader.getCurrentSequenceLength());
    status = POSITION_NOT_FOUND;
  } else {
    // The internal reader is set to the closest bookmark before the
    // position to reach.
    if (position > 0) {
      // If current bookmarked position is not the beginning of the
      // sequence, then goes to the 1-mer preceeding the wanted
      // position.
      if (_reader.getKmerAt(position - 1).empty()) {
        // The wanted position doesn't exist.
        if (sequence.nb_nucleotides == (size_t) -1) {
          sequence.nb_nucleotides = _reader.getCurrentSequenceLength();
        } else {
          assert(sequence.nb_nucleotides == _reader.getCurrentSequenceLength());
        }
        status = POSITION_NOT_FOUND;
      }
    assert(_reader.getState().pos != -1); // TEMP
    }
  }
  return status;

}

DNAFileIndex::Status DNAFileIndex::_setInternalReader(const string &sequence_name, size_t position) {
  ASSERT_VALID_CURRENT_FILE(true);
  DNAFileIndex::Status status = SUCCESS;
  _File &file = _files[_current_file_pos];
  _Sequences &sequences = file.sequences;

  // Check if there is a "loaded" bookmark.
  if (_current_sequence_pos != (size_t) -1) {
    // Yes there is.
    ASSERT_VALID_CURRENT_SEQUENCE(false);
    // Check if the currently loaded bookmark corresponds to the
    // currently processed sequence.
    if (sequences[_current_sequence_pos].name != sequence_name) {
      // Partially "unload" the current bookmark.
      _current_sequence_pos = _getSequencePosition(_current_file_pos, sequence_name);
      _current_bookmark_pos = (size_t) -1;
    }
  } else {
    // If no bookmark is loaded for the current file, then check the
    // already indexed sequences.
    _current_sequence_pos = _getSequencePosition(_current_file_pos, sequence_name);
    assert(_current_bookmark_pos == (size_t) -1);
  }

  if ((_current_sequence_pos == (size_t) -1) && !file.completed) {
    // If the wanted sequence is not yet indexed for this file and the
    // file is not yet completed, then set the internal reader to the
    // last indexed position then parse the file until the wanted
    // sequence is found or the end of file occurs.
    if (file.sequences.empty()) {
      status = _newSequence();
    } else {
      _loadLastBookmark();
    }
    if (status == SUCCESS) {
      bool found = false;
      do {
        assert(_reader);
        _indexEndOfCurrentSequence();
        status = _newSequence();
        if (status == SUCCESS) {
          found = (sequences[_current_sequence_pos].name == sequence_name);
        }
      } while ((status == SUCCESS) && !found);
      if (!found) {
        _current_sequence_pos = (size_t) -1;
        _current_bookmark_pos = (size_t) -1;
        assert(file.completed);
        assert(!_reader);
      }
    }
  }

  // Check if the current sequence was found in the file
  if (_current_sequence_pos != (size_t) -1) {
    status = _setInternalReader(position);
  } else {
    // status should be either SEQUENCE_NOT_FOUND, POSITION_NOT_FOUND,
    // SEQUENCE_NAME_DUPLICATED. or FILE_PARSE_ERROR but neither
    // SUCCESS nor FILE_NOT_FOUND.
    assert(status != SUCCESS);
    assert(status != FILE_NOT_FOUND);
    assert(_current_bookmark_pos == (size_t) -1);
    _reader.reset();
  }

  return status;

}

DNAFileIndex::Status DNAFileIndex::_setInternalReader(const string &filename, const string &sequence_name, size_t position) {

  Status status = SUCCESS;

  // Check if there is a loaded bookmark
  if (_current_file_pos != (size_t) -1) {
    ASSERT_VALID_CURRENT_FILE(false);
    // Check if the currently loaded bookmark corresponds to the given
    // filename.
    if (_reader.getFilename() != filename) {
      // "unload" the current bookmark and close the internal reader.
      _current_file_pos = (size_t) -1;
      _current_sequence_pos = (size_t) -1;
      _current_bookmark_pos = (size_t) -1;
      _reader.close();
    }
  }

  // If no bookmark is loaded for the current file, then the internal
  // reader must be (re)opened with the correct filename.
  if (_current_file_pos == (size_t) -1) {
    assert(_current_sequence_pos == (size_t) -1);
    assert(_current_bookmark_pos == (size_t) -1);
    assert(!_reader);
    _current_file_pos = _getFilePosition(filename);
    if (_current_file_pos != (size_t) -1) {
      assert(_current_file_pos < _files.size());
      _reader.open(filename);
      // The file should exist since it was already processed.
      assert(_reader);
    } else {
      status = _newFile(filename);
    }
  }

  if (_current_file_pos != (size_t) -1) {
    assert(status == SUCCESS);
    status = _setInternalReader(sequence_name, position);
  } else {
    if (_reader) {
      _reader.close();
      assert(status == FILE_PARSE_ERROR);
    } else {
      assert(status == FILE_NOT_FOUND);
    }
  }

  return status;

}

void DNAFileIndex::_setReaderFromInternalReader(DNAFileReader &reader) const {
  if (_current_file_pos == (size_t) -1) {
    assert(_current_sequence_pos == (size_t) -1);
    assert(_current_bookmark_pos == (size_t) -1);
    reader.close();
  } else {
    DNAFileReader::FileState state = _reader.getState();
    state.k = reader.k();
    state.kmer = state.kmer_aux = string(state.k, '?');
    reader.setState(state);
  }
}

DNAFileIndex::DNAFileIndex(const size_t bookmark_distance):
  _files(), _file_locations(), _sequence_file_association(),
  _reader(1),
  _current_file_pos(-1), _current_sequence_pos(-1), _current_bookmark_pos(-1),
  bookmark_distance(bookmark_distance) {
  _reader.check_consistency = false;
  _reader.warn = true;
}

DNAFileIndex::Status DNAFileIndex::sync(const DNAFileReader &reader) {
  DNAFileIndex::Status status = FILE_NOT_FOUND;
  if (reader.getCurrentSequenceID()) {
    const string &desc = reader.getCurrentSequenceDescription();
    const string &name = getNameFromDescription(desc);
    const size_t pos = reader.getCurrentSequenceProcessedNucleotides();
    status = _setInternalReader(reader.getFilename(), name, pos);
  } else {
    if (reader) {
      assert(reader.getCurrentSequenceProcessedNucleotides() == 0);
      assert(reader.getCurrentSequenceDescription().empty());
      assert(reader.getCurrentKmer().empty());
      const _FileLocations::const_iterator it = _file_locations.find(reader.getFilename());
      if (it == _file_locations.end()) {
        status = _newFile(reader.getFilename());
      } else {
        _reader.open(reader.getFilename());
        if (_reader) {
          _current_file_pos = it->second;
          status = SUCCESS;
        } else {
          _current_file_pos = (size_t) -1;
        }
        _current_sequence_pos = (size_t) -1;
        _current_bookmark_pos = (size_t) -1;
      }
    }
  }
  return status;
}

DNAFileIndex::Status DNAFileIndex::set(DNAFileReader &reader, const string &filename, const string &sequence_name, size_t position) {
  Status status = _setInternalReader(filename, sequence_name, position);
  _setReaderFromInternalReader(reader);
  return status;
}

DNAFileIndex::Status DNAFileIndex::indexCurrentSequence(DNAFileReader &reader) {
  Status status = sync(reader);
  if (status == SUCCESS) {
    ASSERT_VALID_CURRENT_FILE(true);
    if (_current_sequence_pos == (size_t) -1) {
      assert(_reader.getCurrentKmer().empty());
      assert(_reader.getCurrentSequenceDescription().empty());
      assert(_reader.getCurrentSequenceID() == 0);
      if (_files[_current_file_pos].sequences.empty()) {
        status = _newSequence();
      } else {
        _current_sequence_pos = 0;
        _loadLastBookmark();
      }
    }
    if (status == SUCCESS) {
      _indexEndOfCurrentSequence();
    } else {
      _reader.gotoSequenceEnd();
    }
    _setReaderFromInternalReader(reader);
  }
  return status;
}

DNAFileIndex::Status DNAFileIndex::indexFile(const string &filename) {
  // Sequence name can't contains (nor start) with spaces, thus trying
  // to set the internal reader to some inexistant sequence name will
  // force whole file indexation.
  Status status = _setInternalReader(filename, " fake name", -1);
  // The SEQUENCE_NOT_FOUND status means file indexation
  // success. Other status are errors and SUCCESS status should never
  // happens.
  assert(status != SUCCESS);
  return ((status == SEQUENCE_NOT_FOUND) ? SUCCESS : status);
}

void DNAFileIndex::clear() {
  _files.clear();
  _file_locations.clear();
  _sequence_file_association.clear();
  _current_file_pos = (size_t) -1;
  _current_sequence_pos = (size_t) -1;
  _current_bookmark_pos = (size_t) -1;
}

size_t DNAFileIndex::getNumberOfFiles() const {
  assert(_files.size() == _file_locations.size());
  return _files.size();
}

list<string> DNAFileIndex::getFilenames() const {
  std::list<string> l;
  for (const auto &f: _files) {
    l.push_back(f.filename);
  }
  return l;
}

size_t DNAFileIndex::getNumberOfSequences(const string &filename) const {
  size_t nb = 0;
  if (filename.empty()) {
    nb = _sequence_file_association.size();
  } else {
    size_t p = _getFilePosition(filename);
    if (p != (size_t) -1) {
      assert(p < _files.size());
      nb = _files[p].sequences.size();
    }
  }
  return nb;
}

list<string> DNAFileIndex::getSequenceNames(const string &filename) const {
  std::list<string> l;
  if (filename.empty()) {
    for (const auto &sf: _sequence_file_association) {
      l.push_back(sf.first);
    }
  } else {
    size_t p = _getFilePosition(filename);
    if (p != (size_t) -1) {
      for (const auto &sf: _files[p].sequences) {
        l.push_back(sf.name);
      }
    }
  }
  return l;
}

string DNAFileIndex::getSequenceName(const string &filename, size_t id) const {
  string name;
  const _Sequence *ptr = ((id > 0) ? _getSequence(filename, id - 1) : NULL);
  if (ptr) {
    name = ptr->name;
  }
  return name;
}

string DNAFileIndex::getSequenceDescription(const string &filename, size_t id) const {
  string desc;
  const _Sequence *ptr = ((id > 0) ? _getSequence(filename, id - 1) : NULL);
  if (ptr) {
    desc = ptr->description;
  }
  return desc;
}

string DNAFileIndex::getSequenceDescription(const string &filename, const string &sequence_name) const {
  string desc;
  const _Sequence *ptr = _getSequence(filename, sequence_name);
  if (ptr) {
    desc = ptr->description;
  }
  return desc;
}

string DNAFileIndex::getSequenceDescription(const string &sequence_name) const {
  return getSequenceDescription(getSequenceFile(sequence_name), sequence_name);
}

size_t DNAFileIndex::getSequenceID(const string &filename, const string &sequence_name) const {
  return _getSequencePosition(_getFilePosition(filename), sequence_name) + 1;
}

size_t DNAFileIndex::getSequenceID(const string &sequence_name) const {
  return getSequenceID(getSequenceFile(sequence_name), sequence_name);
}

string DNAFileIndex::getSequenceFile(const string &sequence_name) const {
  string fname;
  const _SequenceFileAssociation::const_iterator it = _sequence_file_association.find(sequence_name);
  if (it != _sequence_file_association.end()) {
    // Sequence name is associated to filename pointed by *it
    assert(it->first == sequence_name);
    fname = it->second;
    assert(_getFilePosition(fname) != (size_t) -1);
  }
  return fname;
}

size_t DNAFileIndex::getNumberOfNucleotides() const {
  size_t nb = 0;
  for (const auto &f: _files) {
    nb += f.nb_nucleotides;
  }
  return nb;
}

size_t DNAFileIndex::getNumberOfNucleotides(const string &filename) const {
  size_t nb = 0;
  const size_t p = _getFilePosition(filename);
  if (p != (size_t) -1) {
    assert(p < _files.size());
    nb = _files[p].nb_nucleotides;
  }
  return nb;
}

size_t DNAFileIndex::getNumberOfNucleotides(const string &filename,
                                            const string &sequence_name) const {
  size_t nb = 0;

  const _Sequence *ptr = _getSequence(filename, sequence_name);
  if (ptr) {
    // The sequence exists
    const _Bookmarks &bookmarks = ptr->bookmarks;
    if (ptr->nb_nucleotides == (size_t) -1) {
      if (!bookmarks.empty()) {
        nb = (bookmarks.size() - 1) * bookmark_distance;
      }
    } else {
      nb = ptr->nb_nucleotides;
    }
  }
  return nb;
}

DNAFileReader::Format DNAFileIndex::getFileFormat(const string &filename) const {
  DNAFileReader::Format format = DNAFileReader::UNDEFINED_FORMAT;
  size_t p = _getFilePosition(filename);
  if (p != (size_t) -1) {
    assert(p < _files.size());
    format = _files[p].format;
  }
  return format;
}

string DNAFileIndex::getNameFromDescription(const string &description) {
  string name;
  size_t p = 0, n = description.length();
  while ((p < n) && (description[p] <= ' ')) {
    ++p;
  }
  if (p < n) {
    size_t nb = 1;
    while ((p + nb < n) && (description[p + nb] > ' ')) {
      ++nb;
    }
    name = description.substr(p, nb);
  }
  return name;
}

END_KIM_NAMESPACE
