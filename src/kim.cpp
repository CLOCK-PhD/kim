/******************************************************************************
*                                                                             *
*  Copyright © 2023-2024 -- IGH / LIRMM / CNRS / UM                           *
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

#include <kim.h>
#include "optionparser.h"

#include <iostream>
#include <map>

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
  COPYRIGHT_OPT,
  VERSION_OPT,
  QUIET_OPT,
  FORCE_OPT,
  INDEX_DIRECTORY_OPT,
  CREATE_INDEX_OPT,
  KMER_LENGTH_OPT,
  KMER_PREFIX_LENGTH_OPT,
  REFERENCE_OPT,
  VARIANTS_OPT,
  OUTPUT_OPT
};

#define KIM_DEFAULT_INDEX_DIRECTORY     "kim_index"
#define KIM_DEFAULT_KMER_LENGTH         27
#define KIM_DEFAULT_KMER_PREFIX_LENGTH  6

#define _str(x) #x
#define stringify(x) _str(x)

const option::Descriptor usage[] =
  {
   { UNKNOWN_OPT,            0, "" , "",                   Arg::None,
     PACKAGE " version " VERSION " -- " PACKAGE_DESCRIPTION "\n\n"
     "Usages:\n"
     "  kim [options] --create-index\n"
     "  kim [options] <file> [<file> ...]\n"
     "  kim [information option]\n"
     "\n" },
   ////////////////////////////////////////////////////////////////////////
   { UNKNOWN_OPT,            0, "" , "",                   Arg::None,
     "Available information options:" },
   { HELP_OPT,               0, "h", "help",               Arg::None,
     "  -h | --help \tPrint usage and exit." },
   { VERSION_OPT,            0, "v", "version",            Arg::None,
     "  -v | --version \tPrint version and exit." },
   { COPYRIGHT_OPT,          0, "" , "copyright",          Arg::None,
     "       --copyright \tPrint copyright summary and exit." },
   { COPYRIGHT_OPT,          1, "" , "full-copyright",     Arg::None,
     "       --full-copyright \tPrint full copyright and exit." },
   ////////////////////////////////////////////////////////////////////////
   { UNKNOWN_OPT,            0, "" , "",                   Arg::None,
     "\n"
     "Options available for both index creation and querying:" },
   { QUIET_OPT,              0, "q", "quiet",              Arg::None,
     "  -q | --quiet \tDon't produce warning or extra informations on"
     " standard error channel." },
   { FORCE_OPT,              0, "f", "force",              Arg::None,
     "  -f | --force \tForce overwriting existing index directory (on index"
     " creation) or output file (on querying). This is dangerous, you are"
     " advertised." },
   { INDEX_DIRECTORY_OPT,    0, "d", "index-dir",          Arg::Required,
     "  -d | --index-dir <dir> \tDirectory containing the index files"
     " (default: " KIM_DEFAULT_INDEX_DIRECTORY ")." },
   { UNKNOWN_OPT,            0, "" , "",                   Arg::None,
     "\n"
     "Options available only for index creation:" },
   { CREATE_INDEX_OPT,       0, "c", "create-index",       Arg::None,
     "  -c | --create-index \tFlag to run in index creation mode. If this"
     " flag is given, the --reference and --variants options must be"
     " provided." },
   { KMER_LENGTH_OPT,        0, "k", "kmer-length",        Arg::Numeric,
     "  -k | --kmer-length <length> \tLength of the k-mers"
     " (default: " stringify(KIM_DEFAULT_KMER_LENGTH) ")." },
   { KMER_PREFIX_LENGTH_OPT, 0, "p", "kmer-prefix-length", Arg::Numeric,
     "  -p | --kmer-prefix-length <length> \tLength of the k-mers prefix"
     " (default: " stringify(KIM_DEFAULT_KMER_PREFIX_LENGTH) ")." },
   { REFERENCE_OPT,          0, "r", "reference",          Arg::Required,
     "  -r | --reference <file> \tReference sequence from which k-mers are"
     " built. It is possible to use several references. Each given reference"
     " must have exactly one variants file associated to it (see option"
     " --variants). Reference file must be fasta or fastq formatted." },
   { VARIANTS_OPT,           0, "v", "variants",           Arg::Required,
     "  -v | --variants <file> \tVariants to index. There must be exactly"
     " one variant file by reference provided (see option --reference). The"
     " variants file must be VCF formatted (specification 4.3 or above). If"
     " the '##reference' metadata is given and matches one of the reference"
     " file basename, then the given variants file is associated to this"
     " reference file. If some variants file hasn't any '##reference'"
     " metadata or this metadata doesn't correspond to any given reference"
     " file, the kim program will consider that the variants file order"
     " follows the unbound sequence file order." },
   ////////////////////////////////////////////////////////////////////////
   { UNKNOWN_OPT,            0, "" , "",                   Arg::None,
     "\n"
     "Options available only for Index querying:" },
   { OUTPUT_OPT,             0, "o", "output",             Arg::Required,
     "-o | --output <file> \tPath to the filename where results are stored"
     " (instead of being prompted to the standard output)." },
   ////////////////////////////////////////////////////////////////////////
   { UNKNOWN_OPT,            0, "" ,  "",                  Arg::None,
     "\n"
     "The kim program has two major features:\n"
     "- The first one is to create a k-mer index related to some variant of"
     " interest.\n"
     "- The second one is to analyse some given input files (containing raw"
     " biological data from some experiment) and to query an existing index"
     " to exhibit variants that are present in the input data.\n"
     "\nTo create an index, the kim program needs a (set of) reference"
     " genome(s) or transcriptome(s) and the variant tof interest (in VCF"
     " format). An example of command line to create an index is:\n"
     "\n  kim --create-index \\"
     "\n      --kmer-length 27 --prefix-length 6 \\"
     "\n      --reference file1.fasta --reference file2.fasta \\"
     "\n      --variants variants_file.vcf \\"
     "\n      --index-dir /where/to/store/index/"
     "\n"
     "\nTo query an existing index, the kim program needs an index and some"
     " files to analyse. An example of command line to query some index is:\n"
     "\n  kim --index-dir /path/to/my/index/ file1.fastq file2.fastq"
     "\n"
     "\nOf course, it is possible to create the index then query it in a"
     " unique command by combining the necessary options."
     "\nThe biological sequence files (input files to analyse and"
     " references) must be either fasta or fastq formatted."
     "\n" },
     {0, 0, 0, 0, 0, 0}
  };


// Dump the file content (up to the number of given characters) to the
// given stream.
void dump_file(const string &filename, size_t n = size_t(-1), ostream &os = cout) {
  ifstream ifs(filename);
  static const size_t buffer_size = 1024;
  char buffer[buffer_size + 1];
  buffer[buffer_size] = '\0';
  while (n && ifs.read(buffer, buffer_size)) {
    if (n >= buffer_size) {
      n -= buffer_size;
    } else {
      buffer[n] = '\0';
      n = 0;
    }
    os << buffer;
  }
  ifs.close();
  if (n) {
    if (n > (size_t) ifs.gcount()) {
      n = ifs.gcount();
    }
    buffer[n] = '\0';
    os << buffer;
  }
}

void show_copyright(bool full) {
  cout << PACKAGE " version " VERSION " -- " PACKAGE_DESCRIPTION "\n\n";
  string fname_orig = "LICENSE.md";
  string fname = FileReader::findFile(fname_orig);
  if (!fname.empty()) {
    dump_file(fname, 508);
  } else {
    cerr << "File '" << fname_orig << "' not found or not readable."
         << " Please ensure the kim program is correctly installed."
         << endl;
  }
  if (full) {
    cout << endl
         << "You will find below the english version of the licence then the french version."
         << endl
         << endl;
    fname_orig = "LICENSE-en.md";
    fname = FileReader::findFile(fname_orig);
    if (!fname.empty()) {
      dump_file(fname);
    } else {
      cerr << "File '" << fname_orig << "' not found or not readable."
           << " Please ensure the kim program is correctly installed."
           << endl;
    }
    cout << endl;
    fname_orig = "LICENSE-fr.md";
    fname = FileReader::findFile(fname_orig);
    if (!fname.empty()) {
      dump_file(fname);
    } else {
      cerr << "File '" << fname_orig << "' not found or not readable."
           << " Please ensure the kim program is correctly installed."
           << endl;
    }
  }
  cout << endl;

  cout << "The kim program also takes benefit from the BitMagic library:" << endl
       << bm::_copyright<true>::_p << endl
       << "This library is released under the terms of the Apache Licence 2.0:" << endl
       << "http://www.apache.org/licenses/LICENSE-2.0" << endl
       << endl;
}

void dump_results(const map<string, VariantIdentification> &results,
                  const KmerVariantGraph &index,
                  ostream &os) {
  os << "---" << endl
     << "SNPs:" << endl;
  for (map<string, VariantIdentification>::const_iterator it = results.cbegin();
       it != results.cend();
       ++it) {
    os << "  - id: " << it->first << endl
       << "    reads:" << endl;
    list<VariantIdentification::ReadID_type> reads = it->second.getReads();
    for (list<VariantIdentification::ReadID_type>::const_iterator read_it = reads.cbegin();
         read_it != reads.cend();
         ++read_it) {
      os << "      - \"" << *read_it << "\": "
         << it->second.getReadScore(*read_it) * index.getVariantCount(it->first)
         << endl;
    }
  }
}

struct _OptionHandler {
  option::Stats  stats;
  option::Option *options;
  option::Option *buffer;
  option::Parser parse;

  _OptionHandler(int argc, char **argv):
    stats(true, usage, argc, argv),
    options(new option::Option[stats.options_max]),
    buffer(new option::Option[stats.buffer_max]),
    parse(usage, argc, argv, options, buffer) {
  };
  ~_OptionHandler() {
    delete [] options;
    delete [] buffer;
  }
};

int main(int argc, char **argv) {

  argc -= (argc>0);
  argv += (argc>0); // skip program name argv[0] if present
  _OptionHandler opts(argc, argv);

  Settings settings;

  if (opts.parse.error()) {
    return 1;
  }

  if (opts.options[HELP_OPT] || (argc == 0)) {
    option::printUsage(cout, usage);
    return 0;
  }

  if (opts.options[VERSION_OPT]) {
    cout << "kim version " << VERSION << endl;
    return 0;
  }

  if (opts.options[COPYRIGHT_OPT]) {
    show_copyright(opts.options[COPYRIGHT_OPT].type());
    return 0;
  }

  try {

    settings.warn(!opts.options[QUIET_OPT]);

    if (settings.warn()) {
      cerr << "kim version " << VERSION << endl;
    }

    if (opts.options[UNKNOWN_OPT]) {
      Exception e("The kim program doesn't accept the following options:");
      for (option::Option* opt = opts.options[UNKNOWN_OPT]; opt; opt = opt->next()) {
        e << " ";
        if (opt->name[0] != '-') {
          e << "-";
        }
        e << opt->name;
      }
      throw e;
    }

    const char *index_directory = KIM_DEFAULT_INDEX_DIRECTORY;
    if (opts.options[INDEX_DIRECTORY_OPT]) {
      index_directory = opts.options[INDEX_DIRECTORY_OPT].arg;
    }

    if (!opts.options[CREATE_INDEX_OPT]) {
      // If not in index creation mode, then must be in query mode.
      if (opts.parse.nonOptionsCount() == 0) {
        throw Exception("Bad usage");
      }
    }

    // Check if all given input files (if any) are readable
    for (int i = 0; i < opts.parse.nonOptionsCount(); ++i) {
      ifstream ifs(opts.parse.nonOption(i));
      if (ifs) {
        ifs.close();
      } else {
        Exception e;
        e << "Unable to find or open the file '"
          << opts.parse.nonOption(i) << "'.";
        throw e;
      }
    }

    // Ensure that given output filename (if provided) is writable.
    if (opts.options[OUTPUT_OPT]) {
      if (!opts.options[FORCE_OPT]) {
        ifstream ifs(opts.options[OUTPUT_OPT].arg);
        if (ifs) {
          ifs.close();
          Exception e;
          e << "File '" << opts.options[OUTPUT_OPT].arg
            << "' already exists. Please use a new file name or use option '--force'.";
          throw e;
        }
      }
      ofstream ofs(opts.options[OUTPUT_OPT].arg);
      if (ofs) {
        ofs.close();
      } else {
        Exception e;
        e << "Unable to open '" << opts.options[OUTPUT_OPT].arg
          << "' file to print results.";
        throw e;
      }
    }

    if (opts.options[CREATE_INDEX_OPT]) {
      if (opts.options[REFERENCE_OPT].count() != opts.options[VARIANTS_OPT].count()) {
        Exception e;
        e << "The number of reference files provided for index creation ("
          << opts.options[REFERENCE_OPT].count()
          << ") must be the same as the number of variants files provided ("
          << opts.options[VARIANTS_OPT].count()
          << ").";
        throw e;
      }

      size_t kmer_length = KIM_DEFAULT_KMER_LENGTH;
      if (opts.options[KMER_LENGTH_OPT]) {
        kmer_length = stoull(opts.options[KMER_LENGTH_OPT].arg);
      }
      settings.setKmerLength(kmer_length);

      size_t kmer_prefix_length = KIM_DEFAULT_KMER_PREFIX_LENGTH;
      if (opts.options[KMER_PREFIX_LENGTH_OPT]) {
        kmer_prefix_length = stoull(opts.options[KMER_PREFIX_LENGTH_OPT].arg);
      }
      settings.setKmerPrefixLength(kmer_prefix_length);
      settings.freeze();

    }

    // If in query mode...
    if (opts.parse.nonOptionsCount()) {

      // Loading the index.
      KmerVariantGraph kim_index(index_directory, settings);
      kim_index.freeze();
      map<string, VariantIdentification> variants_map;
      // Process each input file
      for (int i = 0; i < opts.parse.nonOptionsCount(); ++i) {
        if (settings.warn()) {
          cerr << "Processing file '" << opts.parse.nonOption(i) << "'" << endl;
        }
        FastqFileReader reader(opts.parse.nonOption(i), settings);
        // Process each k-mer from the current input file
        for (string kmer = reader.getNextKmer(); !kmer.empty(); kmer = reader.getNextKmer()) {
          //cout << "K-mer : " << kmer << endl;
          list<KmerVariantEdgesSubindex::KmerVariantAssociation> assoc_variant = kim_index.search(kmer);
          for (list<KmerVariantEdgesSubindex::KmerVariantAssociation>::const_iterator it = assoc_variant.cbegin(); it != assoc_variant.cend(); ++it) {
            // cerr << "Current variant association concerns variant '" << it->variant_node->first << "' which concerns " << it->variant_node->second << " k-mers" << endl;
            // cerr << "This variant is potentially seen in read " << reader.getCurrentRead()
            //      << " by the k-mer at position " << reader.getCurrentKmerRelativePosition()
            //      << " in the read" << endl;
            VariantIdentification &v_ident = variants_map[it->variant_node.variant];
            v_ident.add(reader.getCurrentRead(), reader.getCurrentKmerRelativePosition());
          }
        }
      }

      // Write the results
      if (opts.options[OUTPUT_OPT]) {
        ofstream ofs(opts.options[OUTPUT_OPT].arg);
        if (ofs) {
          dump_results(variants_map, kim_index, ofs);
          ofs.close();
        } else {
          Exception e;
          e << "Unable to open '"
            << opts.options[OUTPUT_OPT].arg
            << "' file to print results.";
          throw e;
        }
      } else {
        dump_results(variants_map, kim_index, cout);
      }
    }

  } catch (const exception &e) {
    cerr << "The program has encountered the following error:" << endl
         << endl
         << " => " << e.what() << endl
         << endl;
    //    option::printUsage(cerr, usage);
    return 1;
  };

  if (settings.warn()) {
    cerr << "That's All, Folks!!!" << endl;
  }

  return 0;
}
