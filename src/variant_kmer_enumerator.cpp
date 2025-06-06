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

#include "variant_kmer_enumerator.h"

#include "config.h"

#include <cassert>

using namespace std;

BEGIN_KIM_NAMESPACE

size_t VariantKmerEnumerator::_missing = 0;
string VariantKmerEnumerator::_chrom = "";
string VariantKmerEnumerator::_fname = "";
unique_ptr<DNAFileReader> VariantKmerEnumerator::_reader;
unique_ptr<DNAFileIndex> VariantKmerEnumerator::_index;
unique_ptr<Settings> VariantKmerEnumerator::_settings;

void VariantKmerEnumerator::init(const Settings &settings,
                                 const vector<string> &dna_files,
                                 size_t bookmark_size) {
  _settings = make_unique<Settings>(settings);
  _index = make_unique<DNAFileIndex>(bookmark_size);
  //_index->checkConsistency(true); // This is mandatory.
  _index->warn(_settings->warn());
  _reader = make_unique<DNAFileReader>(_settings->k() - 1);
  _reader->warn = false; // DNA file index is already set to emit warnings or not...
  _reader->check_consistency = false; // DNA file index internal reader checks for consistency

  // Indexing DNA files
  if (_settings->warn()) {
    cerr << "* Indexing reference files (this may take a while...)" << endl;
  }
  for (auto const &f: dna_files) {
    if (_settings->warn()) {
      cerr << "  - Indexing file '" << f << "'" << endl;
    }
    DNAFileIndex::Status status;
    try {
      status = _index->indexFile(f);
      if (status != DNAFileIndex::SUCCESS) {
        throw Exception(DNAFileIndex::statusString(status));
      }
    } catch (const FileReaderParseError &) {
      size_t bad_id = _index->getNumberOfSequences(f);
      cerr << "An error occurs while reading file " << bad_id << "." << endl
           << "The index creation may fail due to this error and/or the created index may be errnoeuous."
           << endl;
      if (!_index->checkConsistency()) {
        cerr << "You should try to re-run index creation using the --enable-consistency-checking option."
             << endl;
      }
    }
    if (_settings->warn()) {
      cerr << "    => "
           << _index->getNumberOfSequences(f) << " new sequences found." << endl;
    }
  }
  if (_settings->warn()) {
    cerr << "  - Total number of indexed sequences: "
         << _index->getNumberOfSequences() << endl;
  }

  // Wen need to have a unique ID when it is missing in the variant file.
  _missing = 0;
  _chrom = "";
  _fname = "";
}

VariantKmerEnumerator::VariantKmerEnumerator(const vcfpp::BcfRecord &v):
  _v(v), _id(v.ID() == "." ? "MI" + to_string(++_missing) : v.ID()),
  _sub_ids(), _left_kmer(), _right_kmer(), _pos(-1), _alt(), _cur_alt(0),
  _kmer()
{

  if (v.CHROM() != _chrom) {
    // If the sequence name changes, the DNA file describing the
    // sequence may differs from the previous one...
    _chrom = v.CHROM();
    _fname = _index->getSequenceFile(_chrom);
  }

  if (_fname.empty()) {
    if (_settings->warn()) {
      cerr << "    "
           << "Skipping variant '" << v.ID()
           << "' since sequence " << _chrom
           << " doesn't corresponds to any of the indexed sequences."
           << endl;
    }
  } else {
    // In vcfpp, variant.POS() == variant.Start() + 1 ; POS() is
    // 1-based and Start is 0-based. In DNAFileIndex, positions
    // are 0-based.
    size_t left_kmer_pos, left_kmer_len;
    if ((size_t) v.Start() >= _settings->k()) {
      left_kmer_pos = v.POS() - _settings->k();
      left_kmer_len = _settings->k() - 1;
    } else {
      left_kmer_pos = 0;
      left_kmer_len = v.Start();
    }
    DNAFileIndex::Status status = _index->set(*_reader, _fname, _chrom, left_kmer_pos);
    if (status == DNAFileIndex::SUCCESS) {
      _reader->getKmerAt(left_kmer_pos);
      _left_kmer = _reader->getCurrentKmer().substr(0, left_kmer_len);
      size_t p = _left_kmer.find_last_not_of("ACGT");
      if (p != string::npos) {
        _left_kmer = _left_kmer.substr(p + 1);
      }
      size_t right_kmer_pos = v.End();
      size_t right_kmer_start = 0;
      size_t nb = _reader->getCurrentSequenceLength();
      if (right_kmer_pos + _settings->k() - 1 > nb) {
        right_kmer_start = right_kmer_pos + _settings->k() - nb - 1;
        right_kmer_pos = nb - _settings->k() + 1;
      }
      _right_kmer = _reader->getKmerAt(right_kmer_pos).substr(right_kmer_start);
      p = _right_kmer.find_first_not_of("ACGT");
      if (p != string::npos) {
        _right_kmer = _right_kmer.substr(0, p);
      }

      string alt = v.ALT();
      // Most of the time there is one or two alternative. Thus
      // reserving room for two is be enough and prevent reallocation.
      _alt.reserve(2);
      size_t start = 0;
      do {
        size_t end = alt.find(',', start);
        string s = alt.substr(start, end - start);
        if (s.find_first_not_of(".ACGT") == string::npos) {
          if (!s.empty()) {
            if (s == ".") {
              _alt.push_back("");
            } else {
              _alt.push_back(s);
            }
          }
        } else {
          _alt.push_back("-");
        }
        start = end + 1;
      } while (start);
      _sub_ids.reserve(_alt.size());
      if (_alt.size() > 1) {
        for (size_t i = 0; i < _alt.size(); ++i) {
          if (_alt[i] == "-") {
            _sub_ids.push_back("-");
          } else {
            _sub_ids.push_back(_id + "_" + to_string(i + 1));
          }
        }
      } else {
        _sub_ids.push_back(_id);
      }
      if (!_alt.empty()) {
        _pos = 0;
        _kmer.reserve(_settings->k());
      }
    } else {
      if (_settings->warn()) {
        cerr << "    "
             << "Unable to retrieve the left variant bounding "
             << (_settings->k() - 1) << "-mers"
             << " at position " << left_kmer_pos << endl;
      }
    }
  }
}

bool VariantKmerEnumerator::nextVariantKmer() {
  if (_pos == (size_t) -1) return false;
  _kmer.clear();
  size_t k = _settings->k();
  string s = _alt[_cur_alt];
  size_t nl = _left_kmer.size();
  size_t nm = s.size();
  size_t nr = _right_kmer.size();
  if (s == "-" || ((nl + nm + nr) < (_pos + k))) {
    // there is no more available k-mer for the current variant alt.
    if (++_cur_alt < _alt.size()) {
      // Processing the next variant alt
      _pos = 0;
      return nextVariantKmer();
    } else {
      // There is no more variant alt
      _pos = (size_t) -1;
      return false;
    }
  }

  _kmer.reserve(k);

  if (_pos < nl) {
    // The current k-mer overlaps the left (k-1)-mer prefix
    assert(_pos + 1 < _settings->k());
    _kmer = _left_kmer.substr(_pos);
    k -= nl - _pos;
  }

  // There must have at least one nucleotide to overlap.
  assert(k > 0);
  if (!s.empty()) {
    // This first overlapping nucleotide must belongs to the variant
    // alt part.
    assert(_pos < nl + s.size());
    s = s.substr(_pos > nl ? _pos - nl : 0, k);
    assert(k >= s.size());
    k -= s.size();
    _kmer += s;
  }

  if (k > 0) {
    // current k-mer overlaps the right (k - 1)-mer.
    assert(k <= nr);
    _kmer += _right_kmer.substr(0, k);
    k = 0;
  }
  assert(k == 0);
  assert(_kmer.size() == _settings->k());
  assert(_kmer.find_first_not_of("ACGT") == string::npos);
  ++_pos;

  return true;
}

KmerVariantGraph::Edge VariantKmerEnumerator::getCurrentKmerVariantEdge() const {
  assert(_kmer.size() == _settings->k());
  assert(_pos > 0);
  return { _kmer, _pos - 1, _sub_ids[_cur_alt] };
}

bool VariantKmerEnumerator::getAlleleFrequencies(vector<float> &af, const string &tag) const {
  bool ok = true;
  vcfpp::BcfRecord v = _v;

  try {
    ok = v.getINFO(tag, af);
  } catch (const invalid_argument &) {
    ok = false;
  }
  if (!ok && (tag == "AF")) {
    // Trying with 'AC' and 'AN' fields
    try {
      vector<int> ac;
      int an = 0;
      ok = v.getINFO("AN", an) && (an != 0) && v.getINFO("AC", ac);
      if (ok) {
        af.resize(ac.size());
        for (size_t i = 0; i < ac.size(); ++i) {
          af[i] = ((float) ac[i]) / (float) an;
        }
      }
    } catch (const invalid_argument &e) {
      ok = false;
    }
  }
  assert(!ok || (af.size() == _sub_ids.size()));
  return ok;
}

END_KIM_NAMESPACE
