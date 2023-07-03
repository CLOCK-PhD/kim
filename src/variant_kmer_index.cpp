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

#include "variant_kmer_index.h"
#include "config.h"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <dirent.h>
#include <sys/types.h>

using namespace std;

BEGIN_KIM_NAMESPACE


/////////////// Some stuff to handling parse error ///////////////

class VariantKmerIndexParseError: public exception {

private:
  string s;

public:

  VariantKmerIndexParseError() {
  }

  inline virtual const char *what() const noexcept {
    return s.c_str();
  }

  template <typename T>
  inline VariantKmerIndexParseError &operator<<(const T &t) {
    stringstream ss;
    ss << t;
    s += ss.str();
    return *this ;
  }

};

#define WARNING_MSG(warn, msg)  \
  if (warn) {                   \
    cerr << "Warning:"          \
         << msg << endl;        \
  }                             \
  (void) 0

#define ERROR_MSG(msg)          \
  VariantKmerIndexParseError e; \
  e << msg;                     \
  throw e

//////////////////////////////////////////////////////////////////


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
  if (v == -1) {
    ERROR_MSG("Unable to encode character '" << c << "'");
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

size_t VariantKmerIndex::checkFilenameCorrectness(const string &filename) {
  if (filename.length() != k1) {
    ERROR_MSG("Name of the index file '" << filename << "' should have length " << k1);
  }
  size_t prefix;
  try {
    prefix = encode(filename, k1);
  } catch (...) {
    ERROR_MSG("Name of the index file '" << filename
              << "' should contains only 'A', 'C', 'G', 'T' or 'U' character.\n");
  }
  return prefix;
}

void VariantKmerIndex::parseFile(const string &filename, size_t prefix, bool is_first) {
  ifstream ifs(filename);
  size_t line = 0;
  if (!ifs) {
    ERROR_MSG("Unable to open index file '" << filename << "'");
  }
  while (ifs) {
    string suffix;
    VariantKmerAssociation assoc;
    string chrom, variation, snp_pos, kmer_pos, kmercount;
    // Lecture de la ligne
    // Fichier actuel (à modifier) : 
    // Suffixe, rsid, chrom, variation, snp_pos, kmer_pos, nbr de kmers générés, ingenome 
    ifs >> suffix >> assoc.rs_id >> chrom >> variation >> snp_pos >> assoc.kmer_rank >> kmercount >>assoc.in_genome;
    // test de vérification - on croise les doigts - ça s'affiche !
    //cout << suffix << " " << assoc.rs_id << " " << chrom << " " << variation << " " << snp_pos<< assoc.kmer_rank << " " << kmercount << assoc.in_genome << endl;
    //cout << suffix << " " << assoc.rs_id << " " << assoc.kmer_rank << " " << assoc.in_genome << endl;
    if (ifs) {
      ++line;
      if (is_first) {
        k2 = suffix.length();
        is_first = false;
      } else {
        if (suffix.length() != k2) {
          ERROR_MSG("Badly formatted index file '" << filename << "' (line " << line << ": suffix '" << suffix << "' should have length " << k2 << ").");
        }
      }
      index[prefix].emplace(suffix, assoc);
    }
  }
  ifs.close();
}

VariantKmerIndex::VariantKmerIndex(const char *path, bool warn):
  k(0), k1(0), k2(0),
  index(),
  warn(warn)
{
  DIR *dir = opendir(path);
  if (dir == NULL) {
    ERROR_MSG("Unable to open the given directory index '" << path << "'");
  }

  struct dirent *entry;
  size_t cpt = 0;
  while ((entry = readdir(dir))) {
    string fname(entry->d_name);
    if (!fname.empty() && (fname[0] != '.')) {
      if (!cpt) {
        k1 = fname.length();
        index.resize(1 << (k1 << 1));
      }
      size_t prefix = checkFilenameCorrectness(fname);
      fname = path;
      fname += "/";
      fname += entry->d_name;
      parseFile(fname, prefix, !cpt);
      ++cpt;
    }
  }
  closedir(dir);
  k = k1 + k2;
  if (cpt != (1 << (k1 << 1))) {
    WARNING_MSG(warn, "There is some missing files in the '" << path << "' index directory.");
  }
}

size_t VariantKmerIndex::getKmerLength() const {
  return k;
}

VariantKmerIndex::VariantKmerAssociation getSecond(const VariantKmerIndex::PartialIndex_type::value_type v) {
  return v.second;
}

list<VariantKmerIndex::VariantKmerAssociation> VariantKmerIndex::search(const string &kmer) const {
  list<VariantKmerAssociation> l;
  size_t prefix = encode(kmer, k1);
  const PartialIndex_type &variant_kmer_assoc = index[prefix];
  pair<PartialIndex_type::const_iterator, PartialIndex_type::const_iterator> range = variant_kmer_assoc.equal_range(kmer.substr(k1));
  transform(range.first, range.second, l.begin(), getSecond);
  return l;
}

END_KIM_NAMESPACE
