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

#include <iostream>
#include <unordered_map>
#include "optionparser.h"
#include "file_reader.h"
#include "variant_kmer_index.h"
#include "variant_identification.h"
#include "config.h"

using namespace std;
using namespace kim;

struct Arg: public option::Arg
{
  static void printError(const char* msg1, const option::Option& opt, const char* msg2)
  {
    fprintf(stderr, "ERROR: %s", msg1);
    fwrite(opt.name, opt.namelen, 1, stderr);
    fprintf(stderr, "%s", msg2);
  }
  static option::ArgStatus Unknown(const option::Option& option, bool msg)
  {
    if (msg) printError("Unknown option '", option, "'\n");
    return option::ARG_ILLEGAL;
  }
  static option::ArgStatus Required(const option::Option& option, bool msg)
  {
    if (option.arg != 0)
      return option::ARG_OK;
    if (msg) printError("Option '", option, "' requires an argument\n");
    return option::ARG_ILLEGAL;
  }
  static option::ArgStatus NonEmpty(const option::Option& option, bool msg)
  {
    if (option.arg != 0 && option.arg[0] != 0)
      return option::ARG_OK;
    if (msg) printError("Option '", option, "' requires a non-empty argument\n");
    return option::ARG_ILLEGAL;
  }
  static option::ArgStatus Numeric(const option::Option& option, bool msg)
  {
    char* endptr = 0;
    if (option.arg != 0 && strtol(option.arg, &endptr, 10)){};
    if (endptr != option.arg && *endptr == 0)
      return option::ARG_OK;
    if (msg) printError("Option '", option, "' requires a numeric argument\n");
    return option::ARG_ILLEGAL;
  }
};

enum  optionIndex {
  UNKNOWN_OPT,
  HELP_OPT,
  VERSION_OPT,
  INDEX_DIRECTORY_OPT
};

const option::Descriptor usage[] =
  {
   { UNKNOWN_OPT,         0, "" , "",          option::Arg::None, ("Usage: kim [options] <file> [<file> ...]\n\n"
                                                                   "Options:") },
   { HELP_OPT,            0, "h", "help",      option::Arg::None, "  -h | --help  \tPrint usage and exit." },
   { VERSION_OPT,         0, "v", "version",   option::Arg::None, "  -v | --version  \tPrint version and exit." },
   { INDEX_DIRECTORY_OPT, 0, "d", "index-dir", Arg::Required,     "  -d | --index-dir <dir>  \tDirectory containing the index files." },
   { UNKNOWN_OPT,         0, "" ,  "",         option::Arg::None, ("\nInput files are expected to be fastq formatted\n"
                                                                   "Example:\n"
                                                                   "  kim --index-dir /path/to/my/index/ file1.fastq file2.fastq\n") },
   {0, 0, 0, 0, 0, 0}
  };

int main(int argc, char **argv) {

  argc -= (argc>0);
  argv += (argc>0); // skip program name argv[0] if present
  option::Stats  stats(true, usage, argc, argv);
  option::Option options[stats.options_max],
                 buffer[stats.buffer_max];
  option::Parser parse(usage, argc, argv, options, buffer);

  if (parse.error()) {
    return 1;
  }

  if (options[HELP_OPT] || (argc == 0)) {
    option::printUsage(std::cout, usage);
    return 0;
  }

  if (options[VERSION_OPT]) {
    cout << "kim version " << VERSION << endl;
    return 0;
  }

  for (option::Option* opt = options[UNKNOWN_OPT]; opt; opt = opt->next()) {
    std::cout << "Unknown option: " << opt->name << "\n";
  }

  if (options[UNKNOWN_OPT]) {
    option::printUsage(std::cout, usage);
    return 1;
  }

  const char *index_directory = NULL;
  if (options[INDEX_DIRECTORY_OPT]) {
    index_directory = options[INDEX_DIRECTORY_OPT].arg;
  } else {
    index_directory = "kim_index";
  }

  try {
    cout << "Index directory: '" << index_directory << "'" << endl;
    cout << "Loading index..." << endl;
    VariantKmerIndex kim_index(index_directory);
    cout << "Index loaded" << endl;
    size_t k = kim_index.getKmerLength()+10;
    unordered_map<VariantKmerIndex::VariantID_type, VariantIdentification> variants_map;

    // Process each input file
    for (int i = 0; i < parse.nonOptionsCount(); ++i) {
      cout << "Processing file '" << parse.nonOption(i) << "'" << endl;
      FileReader reader(parse.nonOption(i), k);
      // Process each k-mer from the current input file
      for (std::string kmer = reader.getNextKmer(); !kmer.empty(); kmer = reader.getNextKmer()) {
        list<VariantKmerIndex::VariantKmerAssociation> variant_ids = kim_index.search(kmer);
        for (list<VariantKmerIndex::VariantKmerAssociation>::const_iterator it = variant_ids.cbegin(); it != variant_ids.cend(); ++it) {
          VariantIdentification &v_ident = variants_map[it->rs_id];
          v_ident.add(reader.getCurrentRead(), reader.getCurrentKmerRelativePosition());
        }
      }
    }

    cout << "---" << endl
         << "SNPs:" << endl;
    for (unordered_map<VariantKmerIndex::VariantID_type, VariantIdentification>::const_iterator it = variants_map.cbegin();
         it != variants_map.cend();
         ++it) {
      cout << "  - id: " << it->first << endl
           << "    reads:" << endl;
      list<VariantIdentification::ReadID_type> reads = it->second.getReads();
      for (list<VariantIdentification::ReadID_type>::const_iterator read_it = reads.begin();
           read_it != reads.end();
           ++read_it) {
        cout << "      - " << *read_it << ":" << it->second.getReadScore(*read_it) << endl;
      }
    }
  } catch (const exception &e) {
    cerr << "The program has encountered the following error:" << endl
         << " => " << e.what() << endl;
  };

  cout << "That's All, Folks!!!" << endl;
  return 0;
}
