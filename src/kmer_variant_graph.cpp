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

#include "kmer_variant_graph.h"

#include "config.h"

#include <fstream>
#include <iostream>
#include <dirent.h>
#include <sys/types.h>

#include <unistd.h>
#include <cassert>

using namespace std;

BEGIN_KIM_NAMESPACE


#define WARNING_MSG(msg)        \
  if (_settings.warn()) {       \
    cerr << "Warning:"          \
         << msg << endl;        \
  }                             \
  (void) 0

#define ERROR_MSG(kind, msg)                    \
  do {                                          \
    KmerVariantGraph ## kind error;             \
    error << msg;                               \
    throw error;                                \
  } while (0)

#define PARSE_ERROR_MSG(msg) ERROR_MSG(ParseError, msg)

#define CHECK_FROZEN_STATE(expected_state, mth)                         \
  if (!(expected_state)) {                                              \
    ERROR_MSG(Exception,                                                \
              "Graph must be " << (expected_state ? "frozen" : "unfrozen") \
              << " before calling KmerVariantGraph::" << #mth           \
              << "() method.");                                         \
  }                                                                     \
  (void) 0

//////////////////////////////////////////////////////////////////

/////////////// Some stuff to handle ATGC <=> bits ///////////////

size_t encode(char c) {
  size_t v = ((c == 'A')
              ? 0
              : ((c == 'C')
                 ? 1
                 : ((c == 'G')
                    ? 2
                    : ((c == 'T')
                       ? 3
                       : -1))));
  if (v == size_t(-1)) {
    PARSE_ERROR_MSG("Unable to encode character '" << c << "'");
  }
  return v;
}

size_t encode(const string &kmer, size_t k1) {
  size_t v = 0;
  for (size_t i = 0; i < k1; ++i) {
    v <<= 2; // equiv v *= 4;
    v += encode(kmer[i]);
  }
  return v;
}

char decode(size_t v) {
  switch (v) {
  case 0: return 'A';
  case 1: return 'C';
  case 2: return 'G';
  case 3: return 'T';
  default:
    PARSE_ERROR_MSG("Unable to decode value " << (int) v);
  }
}

string decode(size_t v, size_t k1) {
  string s(k1, '?');
  for (size_t i = k1; i > 0; --i) {
    s[i - 1] = decode(v & 3);
    v >>= 2; // equiv v /= 4;
  }
  return s;
}

//////////////////////////////////////////////////////////////////

///////////// Some stuff to handle index file header /////////////

class HeaderInformations {

public:

  enum Tag {
            NB_KMERS_TAG,
            NB_ASSOCIATIONS_TAG,
            MAX_ASSOCIATIONS_TAG
  };

private:

  static map<string, Tag> _string2tag;

  static map<Tag, string> _tag2string;

  size_t _nb_kmers;

  size_t _nb_associations;

  size_t _max_associations;

public:

  static string toString(Tag t) {
    return _tag2string[t];
  }

  static Tag toTag(const string &s) {
    map<string, Tag>::const_iterator it = _string2tag.find(s);
    if (it == _string2tag.end()) {
      string m = "There is no tag associated to string '";
      m += s;
      m += "'";
      throw runtime_error(m);
    }
    return it->second;
  }

  HeaderInformations(): _nb_kmers(-1), _nb_associations(-1), _max_associations(-1) {
    if (_string2tag.empty()) {
      _string2tag["NbKmers"] = NB_KMERS_TAG;
      _string2tag["NbAssociations"] = NB_ASSOCIATIONS_TAG;
      _string2tag["MaxAssociations"] = MAX_ASSOCIATIONS_TAG;
      _tag2string[NB_KMERS_TAG] = "NbKmers";
      _tag2string[NB_ASSOCIATIONS_TAG] = "NbAssociations";
      _tag2string[MAX_ASSOCIATIONS_TAG] = "MaxAssociations";
    }
  }

  operator bool() const {
    return ((_nb_kmers != size_t(-1))
            && (_nb_associations != size_t(-1))
            && (_max_associations != size_t(-1)));
  }

  size_t operator[](const string &s) const {
    return operator[](toTag(s));
  }

  size_t operator[](Tag t) const {
    switch (t) {
    case NB_KMERS_TAG:
      return _nb_kmers;
    case NB_ASSOCIATIONS_TAG:
      return _nb_associations;
    case MAX_ASSOCIATIONS_TAG:
      return _max_associations;
    }
  }

  size_t &operator[](const string &s) {
    return operator[](toTag(s));
  }

  size_t &operator[](Tag t) {
    switch (t) {
    case NB_KMERS_TAG:
      return _nb_kmers;
    case NB_ASSOCIATIONS_TAG:
      return _nb_associations;
    case MAX_ASSOCIATIONS_TAG:
      return _max_associations;
    }
    return _nb_kmers; // For suppressing compiler warning
  }

};

map<string, HeaderInformations::Tag> HeaderInformations::_string2tag;

map<HeaderInformations::Tag, string> HeaderInformations::_tag2string;

//////////////////////////////////////////////////////////////////

void KmerVariantGraph::_resizeSubindexes(size_t total, size_t estimated_nb_kmers) {
  _kmer_nodes.clear();
  _kmer_nodes.shrink_to_fit();
  _kmer_nodes.reserve(total);
  _edges.clear();
  _edges.shrink_to_fit();
  _edges.reserve(total);
  for (size_t prefix = 0; prefix < total; ++prefix) {
    _kmer_nodes.emplace_back(estimated_nb_kmers);
    _edges.emplace_back(_kmer_nodes.back(), _variant_nodes, estimated_nb_kmers);
  }
}

size_t KmerVariantGraph::_checkFilenameCorrectness(const string &filename) const {
  if (filename.length() != _settings.getKmerPrefixLength()) {
    PARSE_ERROR_MSG("Name of the index file '" << filename
                    << "' should have length " << _settings.getKmerPrefixLength());
  }
  size_t prefix;
  try {
    prefix = encode(filename, _settings.getKmerPrefixLength());
  } catch (...) {
    PARSE_ERROR_MSG("Name of the index file '" << filename
                    << "' should contains only 'A', 'C', 'G', 'T' or 'U' character.\n");
  }
  return prefix;
}

void KmerVariantGraph::_parseFile(const string &filename, size_t prefix) {
  ifstream ifs(filename);
  if (!ifs) {
    PARSE_ERROR_MSG("Unable to open index file '" << filename << "'");
  }

  // First pass, adding new encountered variants, count number of
  // distinct k-mer suffixes and the nuimber of edges that will be
  // added.
  HeaderInformations header_infos;
  size_t line = 0;

  // Parse file header (lines beginning with '#')
  while (ifs && (ifs.peek() == '#')) {
    string tag;
    size_t value;
    ifs >> tag >> value;
    tag = tag.substr(1, tag.length() - 2);
    try {
      header_infos[tag] = value;
      // Ignoring the end of the line.
      while (ifs.get() != '\n');
      ++line;
    } catch (exception &e) {
      WARNING_MSG("Ignoring error: " << e.what());
    }
  }

  if (!header_infos) {
    PARSE_ERROR_MSG("Header is missing or incomplete.");
  }

  // Get the number of k-mers in this file.
  size_t nb_nodes = header_infos[HeaderInformations::NB_KMERS_TAG];
  size_t nb_edges = header_infos[HeaderInformations::NB_ASSOCIATIONS_TAG];

  KmerNodesSubindex &kmer_nodes = _kmer_nodes[prefix];
  KmerVariantEdgesSubindex &edges = _edges[prefix];

  if (!kmer_nodes.empty()) {
    PARSE_ERROR_MSG("There is existing k-mer nodes related to prefix " << decode(prefix, _settings.getKmerPrefixLength()) << "."
                    << " This situation must not occur!");
  }
  if (!edges.empty()) {
    PARSE_ERROR_MSG("There is existing edges related to k-mer having prefix " << decode(prefix, _settings.getKmerPrefixLength()) << "."
                    << " This situation must not occur!");
  }
  // Reserve enough space for the k-mer nodes sub-index.
  kmer_nodes.reserve(nb_nodes);
  // Reserve enouch space for the k-mer edges sub-index
  edges.reserve(nb_edges);
  // Updates the number of total k-mers.
  _nb_kmers += nb_nodes;
  _nb_edges += nb_edges;

  while (ifs) {

    // Try reading the current line
    string suffix, variant;
    uint64_t rank;
    bool in_reference;
    ifs >> suffix >> variant >> rank >> in_reference;

    if (ifs) {
      ++line;

      // Ensure that the k-mer suffix has size k2 (or set k2 on first
      // kmer suffix).
      if (_settings.frozen()) {
        if (suffix.length() != _settings.getKmerSuffixLength()) {
          PARSE_ERROR_MSG("Badly formatted index file '" << filename
                          << "' (line " << line << ": suffix '" << suffix
                          << "' should have length " << _settings.getKmerSuffixLength() << ").");
        }
      } else {
        _settings.setKmerLength(_settings.getKmerPrefixLength() + suffix.length());
        _settings.freeze();
        assert(_settings.valid());
        BoundedSizeString::setMaximalSize(_settings.getKmerSuffixLength());
      }

      try {
        // Create a temporary k-mer node
        KmerNodesSubindex::KmerNode kmer_node = { suffix, in_reference };
        if (!kmer_nodes.sorted()) {
          BoundedSizeString prev_kmer = kmer_nodes.back().suffix;
          PARSE_ERROR_MSG("All k-mers must be lexicographically sorted within index file '" << filename
                          << "' but k-mer " << decode(prefix, _settings.getKmerPrefixLength())
                          << "." << kmer_node.suffix
                          << " at line " << line << " is less than "
                          << decode(prefix, _settings.getKmerPrefixLength()) << "." << prev_kmer
                          << " at line " << (line - 1) << ".");
        }
        edges.add(kmer_node, variant, rank);
      } catch (const KmerVariantGraphParseError &err) {
        PARSE_ERROR_MSG("Badly formatted index file '" << filename
                        << "' at line " << line << ": " << err.what());
      }
    }
  }

  ifs.close();

}

KmerVariantGraph::KmerVariantGraph(Settings &settings, size_t estimated_nb_kmers):
  _settings(settings),
  _nb_kmers(0), _nb_edges(0),
  _kmer_nodes(), _variant_nodes(),
  _edges()
{
  _settings.freeze();
  size_t mem;
  if (estimated_nb_kmers == size_t(-1)) {
    // Half the available memory
    size_t l1 = sysconf(_SC_AVPHYS_PAGES) / 2 * sysconf(_SC_PAGESIZE);
    // 1 GB
    size_t l2 = 1 << 30;
    mem = (l2 < l1 ? l2 : l1);
    estimated_nb_kmers = mem / (_settings.getKmerSuffixLength() + 1 + sizeof(BoundedSizeString));
  }
  size_t total = 1 << (_settings.getKmerPrefixLength() << 1);
  estimated_nb_kmers /= total;
  _resizeSubindexes(total, estimated_nb_kmers);
}

KmerVariantGraph::KmerVariantGraph(const char *path, Settings &settings):
  _settings(settings),
  _nb_kmers(0), _nb_edges(0),
  _kmer_nodes(), _variant_nodes(),
  _edges()
{
  DIR *dir = opendir(path);
  if (dir == NULL) {
    PARSE_ERROR_MSG("Unable to open the given directory index '" << path << "'");
  }

  struct dirent *entry;
  size_t total = 0;
  size_t cpt = 0;
  if (_settings.warn()) {
    cerr << "Index directory: '" << path << "'" << endl
         << "Loading index..." << endl
         << "Step 1: Scanning index for k-mers and variants" << endl;
  }
  _settings.unfreeze();
  _settings.setKmerLength(size_t(-1));
  while ((entry = readdir(dir))) {
    string fname(entry->d_name);
    if (!fname.empty() && (fname[0] != '.')) {
      if (!cpt) {
        _settings.setKmerPrefixLength(fname.length());
        total = 1 << (_settings.getKmerPrefixLength() << 1);
        _resizeSubindexes(total, 0);
      }
      size_t prefix = _checkFilenameCorrectness(fname);
      fname = path;
      if (fname.back() != '/') {
        fname += "/";
      }
      fname += entry->d_name;
      if (_settings.warn()) {
        cerr << "\033[sProcessing file "
             << (cpt + 1) << " / " << total
             << ": '" << fname << "'";
      }
      _parseFile(fname, prefix);
      ++cpt;
      if (_settings.warn()) {
        cerr << " (" << (cpt * 10000 / total) / 100. << "%)\033[K\033[u";
      }
      // if (prefix == 284) {
      //   break;
      // }
    }
  }
  if (_settings.warn()) {
    cerr << "\033[K"
         << "Number of indexed variants: " << getNbVariants() << endl
         << "Number of indexed k-mers: " << getNbKmers() << endl
         << "Number of indexed edges: " << getNbKmerVariantEdges() << endl
         << "Index loaded" << endl;
  }

  if (cpt != total) {
    WARNING_MSG("There is some missing files in the '" << path << "' index directory.");
  }

  closedir(dir);
}

void KmerVariantGraph::freeze() {
  for (size_t prefix = 0; prefix < _kmer_nodes.size(); ++prefix) {
    KmerVariantEdgesSubindex &edges = _edges[prefix];
    edges.compact();
    edges.freeze();
  }
}

list<KmerVariantEdgesSubindex::KmerVariantAssociation> KmerVariantGraph::search(const std::string &kmer) const {
  CHECK_FROZEN_STATE(frozen(), search);
  size_t prefix = encode(kmer, _settings.getKmerPrefixLength());
  const KmerVariantEdgesSubindex &edges = _edges[prefix];

  BoundedSizeString suffix = kmer.substr(_settings.getKmerPrefixLength());
  return edges.getKmerVariantAssociation(suffix);
}

KmerVariantGraph &KmerVariantGraph::add(const string &variant, const string &kmer, size_t rank) {
  CHECK_FROZEN_STATE(!frozen(), add);
  size_t prefix = encode(kmer, _settings.getKmerPrefixLength());
  BoundedSizeString suffix = kmer.substr(_settings.getKmerPrefixLength());
  KmerVariantEdgesSubindex &edges = _edges[prefix];
  KmerNodesSubindex::KmerNode kmer_node = { suffix, false };
  edges.add(kmer_node, variant, rank);
  return *this;
}

bool KmerVariantGraph::addKmerNode(const string &kmer, bool in_reference) {
  CHECK_FROZEN_STATE(!frozen(), setInReferenceKmer);
  size_t prefix = encode(kmer, _settings.getKmerPrefixLength());
  KmerNodesSubindex::KmerNode kmer_node = { kmer.substr(_settings.getKmerPrefixLength()), in_reference};
  KmerNodesSubindex &kmer_nodes = _kmer_nodes[prefix];
  return kmer_nodes.add(kmer_node);
}

END_KIM_NAMESPACE
