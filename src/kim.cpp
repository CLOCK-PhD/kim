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

#include <kim.h>

#include <optionparser.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <vcfpp.h>
#pragma GCC diagnostic pop

#include <cmath>
#include <filesystem>
#include <iostream>
#include <map>
#include <regex>

#ifdef _OPENMP
#  include <omp.h>
#endif

using namespace std;
using namespace kim;
namespace fs = filesystem;

struct Arg: public option::Arg {

  static void printError(const char* msg1, const option::Option& opt, const char* msg2);

  static option::ArgStatus Unknown(const option::Option& option, bool msg);

  static option::ArgStatus Required(const option::Option& option, bool msg);

  static option::ArgStatus NonEmpty(const option::Option& option, bool msg);

  static option::ArgStatus Numeric(const option::Option& option, bool msg);

  static option::ArgStatus Percentage(const option::Option& option, bool msg);

};

void Arg::printError(const char* msg1, const option::Option& opt, const char* msg2) {
  fprintf(stderr, "ERROR: %s", msg1);
  fwrite(opt.name, opt.namelen, 1, stderr);
  fprintf(stderr, "%s", msg2);
}

option::ArgStatus Arg::Unknown(const option::Option& option, bool msg) {
  if (msg) printError("Unknown option '", option, "'\n");
  return option::ARG_ILLEGAL;
}

option::ArgStatus Arg::Required(const option::Option& option, bool msg) {
  if (option.arg != NULL) return option::ARG_OK;
  if (msg) printError("Option '", option, "' requires an argument\n");
  return option::ARG_ILLEGAL;
}

option::ArgStatus Arg::NonEmpty(const option::Option& option, bool msg) {
  if ((option.arg != NULL) && (option.arg[0] != '\0')) return option::ARG_OK;
  if (msg) printError("Option '", option, "' requires a non-empty argument\n");
  return option::ARG_ILLEGAL;
}

option::ArgStatus Arg::Numeric(const option::Option& option, bool msg) {
  char* endptr = NULL;
  if (option.arg != NULL) strtol(option.arg, &endptr, 10);
  if ((endptr != option.arg) && (*endptr == '\0')) return option::ARG_OK;
  if (msg) printError("Option '", option, "' requires a numeric argument\n");
  return option::ARG_ILLEGAL;
}

option::ArgStatus Arg::Percentage(const option::Option& option, bool msg) {
  char* endptr = NULL;
  double v;
  if (option.arg != NULL) {
    v = strtod(option.arg, &endptr);
    if (endptr != option.arg) {
      if (*endptr == '%') {
        ++endptr;
        v *= 0.01;
      }
      if ((*endptr == '\0') && (v >= 0) && (v <= 1)) return option::ARG_OK;
    }
  }
  if (msg) printError("Option '", option, "' requires a value in [0;1] (or 0% and 100%)\n");
  return option::ARG_ILLEGAL;
}

// A record of the current read description and ID for read analysis callbacks
struct _currentReadInfos {
  string description;
  size_t id;
};

class KimProgram {

private:

  static const option::Descriptor _OPTION_EMPTY_LINE;
  static const option::Descriptor _OPTION_PREAMBLE;
  static const option::Descriptor _OPTION_INFORMATIONS_HEADER;
  static const option::Descriptor _OPTION_INFORMATIONS_HELP;
  static const option::Descriptor _OPTION_INFORMATIONS_VERSION;
  static const option::Descriptor _OPTION_INFORMATIONS_COPYRIGHT;
  static const option::Descriptor _OPTION_COMMON_OPTIONS_HEADER;
  static const option::Descriptor _OPTION_COMMON_OPTIONS_QUIET;
  static const option::Descriptor _OPTION_COMMON_OPTIONS_FORCE;
  static const option::Descriptor _OPTION_COMMON_OPTIONS_INDEX_DIRECTORY;
  static const option::Descriptor _OPTION_COMMON_OPTIONS_CHECK_CONSISTENCY;
  static const option::Descriptor _OPTION_INDEX_CREATION_OPTIONS_HEADER;
  static const option::Descriptor _OPTION_INDEX_CREATION_OPTIONS_CREATE_INDEX;
  static const option::Descriptor _OPTION_INDEX_CREATION_OPTIONS_KMER_LENGTH;
  static const option::Descriptor _OPTION_INDEX_CREATION_OPTIONS_KMER_PREFIX_LENGTH;
  static const option::Descriptor _OPTION_INDEX_CREATION_OPTIONS_KMER_FILTER;
  static const option::Descriptor _OPTION_INDEX_CREATION_OPTIONS_REFERENCE;
  static const option::Descriptor _OPTION_INDEX_CREATION_OPTIONS_VARIANTS;
  static const option::Descriptor _OPTION_INDEX_CREATION_OPTIONS_VARIANT_FILTER;
  static const option::Descriptor _OPTION_INDEX_CREATION_OPTIONS_ALLELE_FREQUENCY_TAG;
  static const option::Descriptor _OPTION_QUERY_OPTIONS_HEADER;
  static const option::Descriptor _OPTION_QUERY_OPTIONS_OUTPUT_DIR;
  static const option::Descriptor _OPTION_QUERY_OPTIONS_ALPHA;
  static const option::Descriptor _OPTION_QUERY_OPTIONS_THRESHOLD;
  static const option::Descriptor _OPTION_QUERY_OPTIONS_MODE;
  static const option::Descriptor _OPTION_QUERY_OPTIONS_DUMP_READ_ANALYSIS;
  static const option::Descriptor _OPTION_FOOTER_ABOUT;
  static const option::Descriptor _OPTION_FOOTER_INDEX_CREATION;
  static const option::Descriptor _OPTION_FOOTER_QUERY;
  static const option::Descriptor _OPTION_END;

  static const option::Descriptor _short_usage[];
  static const option::Descriptor _create_index_usage[];
  static const option::Descriptor _query_usage[];
  static const option::Descriptor _full_usage[];

  struct _OptionHandler {
    option::Stats  stats;
    option::Option *options;
    option::Option *buffer;
    option::Parser parse;

    _OptionHandler(int argc, char **argv);
    ~_OptionHandler();
  };

  enum  _OptionIndex {
    // Information options
    UNKNOWN_OPT,
    HELP_OPT,
    COPYRIGHT_OPT,
    VERSION_OPT,
    // Common options
    QUIET_OPT,
    FORCE_OPT,
    INDEX_DIRECTORY_OPT,
    CONSISTENCY_CHECKING_OPT,
    // Index creation options
    CREATE_INDEX_OPT,
    VARIANT_FILTER_OPT,
    KMER_LENGTH_OPT,
    KMER_PREFIX_LENGTH_OPT,
    KMER_FILTER_OPT,
    REFERENCE_OPT,
    VARIANTS_OPT,
    ALLELE_FREQUENCY_TAG_OPT,
    // Index query options
    OUTPUT_OPT,
    ALPHA_OPT,
    THRESHOLD_OPT,
    MODE_OPT,
    DUMP_READ_ANALYSIS_OPT
  };

  enum _RunningMode {
    UNDEFINED_MODE,
    QUERY_MODE,
    INDEX_CREATION_MODE
  };

  string _progname;
  Settings _settings;
  bool _dump_read_analysis_opt;
  _RunningMode _running_mode;
  vector<pair<string, regex>>  _kmer_filters;
  vector<string> _variant_filters;
  vector<string> _dna_files;
  vector<fs::path> _variants_files;
  map<string, fstream> _vaf_files;
  fs::path _tmpdir;
  fs::path _output_dir;

  KimProgram(int argv, char **argc);
  ~KimProgram();

  void _processOptions(_OptionHandler &_opts);
  bool _variantOfInterest(vcfpp::BcfRecord &v) const;
  void _runQuery();
  void _createIndex();
  fs::path _buildResultFile(const fs::path &f) const;

  // Creates a temporary directory (to store computed allele
  // frequencies on the fly).
  void _createTemporaryDirectory();

  // Creates temporary VAF files (open output streams)
  void _createTemporaryVafFiles(KmerVariantGraph &index);

  // Dump the allele frequencies of the given variant for all
  // populations tags.
  //
  // Since frequencies are in the range [0;1] instead of using 32 bits
  // float encoding, we choose to encode this range using a 16 bits
  // integer. The first bit is used to denote alternative alleles (the
  // first alternative allele is thus in range [0; 32767) and the
  // others are in [32768; 65535), so the correct allele frequencies
  // for those one are computed by substracting 32768). That way, a
  // final value of 0 denotes a frequency of 0 and a final value of
  // 32766 denotes a frequency of 1. A final value of i denotes a
  // frequency of i/32766, leading to a precision of 3e-5. The special
  // value 32767 is reserved to denode the value which doesn't exist.
  void _dump2temporaryVafFiles(VariantAlleleFrequencies &vaf);

  friend void _readAnalyzerDumpCallback(const KimProgram &prog, ReadAnalyzer &r, _currentReadInfos &infos, ostream &os);
  void _dumpResultsHeader(const string &input_fname, ostream &os) const;
  void _dumpResultsRead(const ReadAnalyzer::VariantKmerRates &res, const _currentReadInfos &infos, ostream &os) const;
  void _dumpResultsVariants(const ReadAnalyzer::VariantScores &results, ostream &os) const;
  static void _dumpResultsKimScores(const map<string, double> &scores, ostream &os);
  static void _dumpResultsFooter(Monitor &monitor, ostream &os);

  // Closes and deletes temporary VAF files
  void _deleteTemporaryDirectory();

  // Saves temporary VAF files to final index directory
  void _saveTemporaryVafFiles();

  // Load VAF files (open input streams) from index directory
  void _loadIndexVafFiles(const KmerVariantGraph &index);

  // Closes opened VAF files (either temporary or from index)
  void _closeVafFiles();

  // Computes kim scores.
  map<string, double> _computeKimScores(const KmerVariantGraph &index, const ReadAnalyzer::VariantScores &result);

  // Check whether the given string is an URL.
  static bool _isURL(const string &f);

  // Check if the given file is readable.
  static bool _checkIfFileIsReadable(const fs::path &f);

  // Ensures that the given file is readable.
  static void _assertFileIsReadable(const fs::path &f);

  // Dump the file content (up to the number of given characters) to
  // the given stream.
  static void _dumpFile(const fs::path &filename, size_t n = size_t(-1), ostream &os = cout);

  // Display the kim program [summary or full content] copyright
  static void _showCopyright(bool full);

  // Compute the common prefix of the given path
  static fs::path _commonPathPrefix(const vector<string> &files);


public:

  KimProgram(const KimProgram&) = delete;

  KimProgram &operator=(const KimProgram &) = delete;

  static KimProgram &App(int argc, char **argv) {
    static KimProgram _mainApp(argc, argv);
    return _mainApp;
  }

  int run();

};

#define KIM_DEFAULT_INDEX_DIRECTORY     "kim_index"
#define KIM_DEFAULT_RESULT_EXTENSION    ".kyr" // Kim Yaml Result
#define KIM_DEFAULT_KMER_LENGTH         27
#define KIM_DEFAULT_KMER_PREFIX_LENGTH  6
#define KIM_DEFAULT_ALPHA               "1%"
#define KIM_DEFAULT_THRESHOLD           "0%"
#define KIM_DEFAULT_MODE                "weak"
static const string VARIANT_ALLELE_FREQUENCIES_BASENAME = "0_vaf";
string _buildVafFilename(const string &tag) {
  return VARIANT_ALLELE_FREQUENCIES_BASENAME + "_" + tag + ".bin";
}

#define _str(x) #x
#define stringify(x) _str(x)

const option::Descriptor KimProgram::_OPTION_EMPTY_LINE = {
  UNKNOWN_OPT, 0, "" , "", Arg::None, "\n"
};

const option::Descriptor KimProgram::_OPTION_PREAMBLE = {
  UNKNOWN_OPT, 0, "" , "", Arg::None,
  PACKAGE " version " VERSION " -- " PACKAGE_SHORT_DESCRIPTION "\n\n"
  "Usages:\n"
  "  kim [options] --create-index\n"
  "  kim [options] <file> [<file> ...]\n"
  "  kim [information option]\n"
};

const option::Descriptor KimProgram::_OPTION_INFORMATIONS_HEADER = {
  UNKNOWN_OPT, 0, "" , "", Arg::None,
  "Available information options:"
};

const option::Descriptor KimProgram::_OPTION_INFORMATIONS_HELP = {
  HELP_OPT, 0, "h", "help", Arg::Optional,
  "  -h | --help[=<short*|create-index|query|full> \t"
  "Print usage and exit. This flag can be given an option to display "
  "additional details on index creation or on query. Default is to "
  "display a short help message."
};

const option::Descriptor KimProgram::_OPTION_INFORMATIONS_VERSION = {
  VERSION_OPT, 0, "v", "version", Arg::None,
  "  -v | --version \t"
  "Print version and exit."
};

const option::Descriptor KimProgram::_OPTION_INFORMATIONS_COPYRIGHT = {
  COPYRIGHT_OPT, 0, "" , "copyright", Arg::Optional,
  "       --copyright[=<short*|full>] \t"
  "Print copyright summary (short or default) or complete (full) and exit."
};


const option::Descriptor KimProgram::_OPTION_COMMON_OPTIONS_HEADER = {
  UNKNOWN_OPT, 0, "" , "", Arg::None,
  "Options available for both index creation and querying:"
};

const option::Descriptor KimProgram::_OPTION_COMMON_OPTIONS_QUIET = {
  QUIET_OPT, 0, "q", "quiet", Arg::None,
  "  -q | --quiet \t"
  "Don't produce warning or extra informations on standard error channel."
};

const option::Descriptor KimProgram::_OPTION_COMMON_OPTIONS_FORCE = {
  FORCE_OPT, 0, "f", "force", Arg::None,
  "  -f | --force \t"
  "Force overwriting existing index directory. "
  "This is dangerous, you are advertised."
};

const option::Descriptor KimProgram::_OPTION_COMMON_OPTIONS_INDEX_DIRECTORY = {
  INDEX_DIRECTORY_OPT, 0, "d", "index-dir", Arg::Required,
  "  -d | --index-dir <dir> \t"
  "Directory containing the index files  (default: " KIM_DEFAULT_INDEX_DIRECTORY ")."
};

const option::Descriptor KimProgram::_OPTION_COMMON_OPTIONS_CHECK_CONSISTENCY = {
  CONSISTENCY_CHECKING_OPT, 0, "", "check-consistency", Arg::None,
  "       --check-consistency \t"
  "Enable consistency checking while parsing DNA sequence files. "
  "This significantly increases the program running time but ensures that "
  "your file are correctly written (useful when something wrong occurs, "
  "disabled by default)"
};

const option::Descriptor KimProgram::_OPTION_INDEX_CREATION_OPTIONS_HEADER = {
  UNKNOWN_OPT, 0, "" , "", Arg::None,
  "Options available only for index creation:"
};

const option::Descriptor KimProgram::_OPTION_INDEX_CREATION_OPTIONS_CREATE_INDEX = {
  CREATE_INDEX_OPT, 0, "c", "create-index", Arg::None,
  "  -c | --create-index \t"
  "Flag to run in index creation mode. If this flag is given, the "
  "--reference and --variants options must be provided."
};

const option::Descriptor KimProgram::_OPTION_INDEX_CREATION_OPTIONS_KMER_LENGTH = {
  KMER_LENGTH_OPT, 0, "k", "kmer-length", Arg::Numeric,
  "  -k | --kmer-length <length> \t"
  "Length of the k-mers (default: " stringify(KIM_DEFAULT_KMER_LENGTH) ")."
};

const option::Descriptor KimProgram::_OPTION_INDEX_CREATION_OPTIONS_KMER_PREFIX_LENGTH = {
  KMER_PREFIX_LENGTH_OPT, 0, "p", "kmer-prefix-length", Arg::Numeric,
  "  -p | --kmer-prefix-length <length> \t"
  "Length of the k-mers prefix (default: " stringify(KIM_DEFAULT_KMER_PREFIX_LENGTH) ")."
};

const option::Descriptor KimProgram::_OPTION_INDEX_CREATION_OPTIONS_KMER_FILTER = {
  KMER_FILTER_OPT, 0, "E", "exclude-kmer", Arg::Required,
  "  -E | --exclude-kmer <expr> \t"
  "Filter out k-mers matching the given expression. When this option is "
  "provided multiple time, any k-mer that matches some of the given "
  "expression is filtered out from the index. For example, to ignore "
  "poly-A and poly-T k-mers, you can provide either the expression "
  "'^(A+)|(T+)$' or provide the two separate expressions '^A+$' and "
  "'^T+$'. The former is less efficient."
};

const option::Descriptor KimProgram::_OPTION_INDEX_CREATION_OPTIONS_REFERENCE = {
  REFERENCE_OPT, 0, "R", "reference", Arg::Required,
  "  -r | --reference <file> \t"
  "Reference sequence from which k-mers are built. It is possible to use "
  "several references and at least one reference file must be provided "
  "(see option --variants). The given reference file must be fasta or "
  "fastq formatted."
};

const option::Descriptor KimProgram::_OPTION_INDEX_CREATION_OPTIONS_VARIANTS = {
  VARIANTS_OPT, 0, "V", "variants", Arg::Required,
  "  -V | --variants <file> \t"
  "Variants to index. The sequence name of each variant must be defined "
  "in exactly one of the provided reference files (see option "
  "--reference). The variants file must be VCF of BCF formatted "
  "(compressed with gzip or not). You also can provide URL instead of "
  "some local file name."
};

const option::Descriptor KimProgram::_OPTION_INDEX_CREATION_OPTIONS_VARIANT_FILTER = {
  VARIANT_FILTER_OPT, 0, "F", "variant-filter", Arg::Required,
  "  -F | --variant-filter <expr> \t"
  "Filter variants according to the given expression. When this option is "
  "provided multiple time, each expression must be satisfyed for the "
  "variant to be added into the index. See filter syntax explanation "
  "below for additional details."
};

const option::Descriptor KimProgram::_OPTION_INDEX_CREATION_OPTIONS_ALLELE_FREQUENCY_TAG = {
  ALLELE_FREQUENCY_TAG_OPT, 0, "T", "allele-frequency-tag", Arg::Required,
  "  -T | --allele-frequency-tag> \t"
  "Set an extra tag to use for variant allele frequencies. During index "
  "creation, the default variant allele frequencies comes from the 'AF' "
  "info tag of, if this field is missing, they are computed from the 'AC' "
  "and 'AN' info tags if available (see VCF specification). Using this "
  "option allows you to provide an additional allele frequency tag "
  "(providing this option several time allows you to add as many new "
  "tags)."
};

const option::Descriptor KimProgram::_OPTION_QUERY_OPTIONS_HEADER = {
  UNKNOWN_OPT, 0, "" , "", Arg::None,
  "Options available only for Index querying:"
};

const option::Descriptor KimProgram::_OPTION_QUERY_OPTIONS_OUTPUT_DIR = {
  OUTPUT_OPT, 0, "o", "output-dir", Arg::Required,
  "  -o | --output-dir <dir> \t"
  "Path to the directory where results are stored (instead of being "
  "prompted to the standard output). One result file is created by input "
  "file with the '" KIM_DEFAULT_RESULT_EXTENSION "' extension. appended."
};

const option::Descriptor KimProgram::_OPTION_QUERY_OPTIONS_ALPHA = {
  ALPHA_OPT, 0, "a", "alpha", Arg::Percentage,
  "  -a | --alpha <value> \t"
  "Type I error (significance) of variant analysis. This is the "
  "probability to reject a variant which really is in the analyzed "
  "data. (default:" KIM_DEFAULT_ALPHA ")."
};

const option::Descriptor KimProgram::_OPTION_QUERY_OPTIONS_THRESHOLD = {
  THRESHOLD_OPT, 0, "t", "threshold", Arg::Percentage,
  "  -t | --threshold <value> \t"
  "Minimum k-mer rate threshold to consider a variant being present in "
  "some read. The threshold must be in the range [0; 1] or expressed as a "
  "percentage (default:" KIM_DEFAULT_THRESHOLD ")."
};

const option::Descriptor KimProgram::_OPTION_QUERY_OPTIONS_MODE = {
  MODE_OPT, 0, "m", "mode", Arg::Required,
  "  -m | --mode <weak*|strict> \t"
  "Set the query mode. A weak mode means that this analyzer uses all "
  "k-mers to [try to] detect variants whereas a strict mode means that "
  "this analyzer uses only k-mers that are not in the reference (default:"
  KIM_DEFAULT_MODE ")."
};

const option::Descriptor KimProgram::_OPTION_QUERY_OPTIONS_DUMP_READ_ANALYSIS = {
  DUMP_READ_ANALYSIS_OPT, 0, "R", "read-analysis", Arg::None,
  "  -R | --read-analysis \t"
  "Dump per read analysis. This requires two passes over the reads (up to "
  "twice the time required without this option) and may lead to very huge "
  "result files. Unless you really need per read analysis, you shouldn't "
  "run queries with this option."
};

const option::Descriptor KimProgram::_OPTION_FOOTER_ABOUT = {
  UNKNOWN_OPT, 0, "" , "", Arg::None,
  "ABOUT KIM\n"
  "=========\n"
  "\n"
  PACKAGE_FULL_DESCRIPTION "\n"
  "The kim program has two major features:\n"
  "- The first one is to create a k-mer index related to some variant of"
  " interest.\n"
  "- The second one is to analyse some given input files (containing raw"
  " biological data from some experiment) and to query an existing index"
  " to exhibit variants that are present in the input data."
};

const option::Descriptor KimProgram::_OPTION_FOOTER_INDEX_CREATION = {
  UNKNOWN_OPT, 0, "" , "", Arg::None,
  "INDEX CREATION\n"
  "==============\n"
  "\n"
  "To create an index, the kim program needs a (set of) reference"
  " genome(s) or transcriptome(s) and the variant tof interest (in VCF"
  " format). An example of command line to create an index is:\n"
  "\n  kim --create-index \\"
  "\n      --kmer-length 27 --kmer-prefix-length 6 \\"
  "\n      --reference file1.fasta --reference file2.fasta \\"
  "\n      --variants variants_file.vcf \\"
  "\n      --index-dir /where/to/store/index/"
  "\n"
  "The biological sequence files must be either fasta or fastq formatted.\n"
  "\n"
  "Variant Filtering\n"
  "-----------------\n"
  "\n"
  "The filter can be a simple expression of the form:\n"
  "   '<attribute> <comparator> <value>'\n"
  "or a more complex expression using logical operators.\n\n"
  "Accepted (case insensive) operators are:\n"
  "  'not', '!' (shortcut for 'not'),\n"
  "  'and', '&&' (shortcut for 'and'),\n"
  "  'or' and '||' (shortcut for 'or').\n"
  "\n"
  "The conjunction ('and') has priority over the disjunction ('or') and the"
  " negation ('not') has priority over the conjunction. It is possible to use"
  " parenthesis to enclose an expression in order to give it a higher priority."
  " For example '<expr1> or not <expr2> and <expr3>' evaluates exactly as"
  " '(<expr1> or (not(<expr2>) and <expr3>))'.\n"
  "\n"
  "The possible <attribute>s can be either VCF fields or some variant status:\n"
  "- Handled (case insensive) VCF fields are:\n"
  "  'CHROM', 'ID', 'POS', 'FILTER' and 'QUAL'.\n"
  "- It can also be some (CASE SENSIVE) key of the 'INFO' field.\n"
  "  The syntax of the <attribute> is then 'INFO:KEY'.\n"
  "- The following (case insensive) variant properties are also accepted:\n"
  "  'Ploidy', 'SNP', 'MSNP', 'SV' and 'Indel'.\n"
  "\n"
  "The available <comparator>s depends on the type of the attribute.\n"
  "\n"
  "If the attribute has a boolean value ('SNP', 'mSNP', 'SV', 'Indel' and some"
  " info tags like 'INFO:1000G' of 'INFO:R3' for example), then only '=' and"
  " '!=' (or '<>') are accepted.\n"
  "\n"
  "If the attribute has a string value ('CHROM', 'ID', 'FILTER' and some"
  " info tags like 'info:VC' for example), then only '=', '~', '!=' (or '<>')"
  " are accepted. The '~' is a regex comparison whereas the other operators"
  " expect exact (in)equality.\n"
  "\n"
  "If the attribute has a numeric value ('POS', 'QUAL', 'Ploidy' and some"
  " INFO tags like 'info:RS' for example), then available operators are"
  " '=', '!=' (or '<>'), '<', '<=', '>', '>='.\n"
  "\n"
  "The compared <value> must be of the correct type according to the"
  "attribute. For boolean values, only 'true', '1' (shortcut for 'true'),"
  " 'false', '0' (shortcut for 'false') are accepted (case insensive).\n"
  "\n"
  "String values must be enclosed by either single quotes (in such case,"
  " inner single quotes must be escaped uising '\\') or by double quotes"
  " (in such case, inner double quotes must be escaped using '\\').\n"
  "\n"
  "Numeric values are either integer like -3, +54, 42, ... or reals"
  " (possibly using a scientific notation) like -3., +5.12, .18, 1.2e-3,"
  " ...\n"
  "\n"
  "Let us illustate filters using some examples:\n"
  "\n"
  "- Filter 'Chrom~\".*\"' accepts all variants (since any name"
  " is accepted for the 'CHROM' attribute.\n"
  "\n"
  "- Filter 'Chrom=\"B\"' accepts only variants on the sequence named 'B'.\n"
  "\n"
  "- Filter 'Pos>123' accepts any variants located after position 123"
  " (starting from 1) in any sequence.\n"
  "\n"
  "- Filter 'SNP=1' accepts only variants that are SNPs (multi-allelic SNP\n"
  " are not accepted).\n"
  "\n"
  "- Filter 'SNP=1 or MSNP=true' accepts only variants that are SNPs,"
  " including multi-allelic SNP.\n"
  "\n"
  "- Filter 'info:RS>=35' accepts variants having the key 'RS' in the info\n"
  "field with a value greater or equalt to 35."
  "\n"
  "- Filter 'info:VC=\"INS\"' accepts variants having the key 'VC' (variant"
  " category) in the info field with an exact value of 'INS' (insertion).\n"
  "\n"
  "- Filter 'Chrom=\"D\" and (pos < 50 or pos > 100)' accepts only variants"
  " located on sequence 'D' before position 50 or after position 100.\n"
  "\n"
  "- Filter 'not (Chrom~\"[A-C]\" or pos >= 50 and pos <= 100)' accepts"
  " variants that are located in any sequence except sequences 'A', 'B' and 'C'."
  " It also accepts any variants located between positions 50 and 100 on any"
  " sequence. The 'and' operator has the priority over the 'or' operator.\n"
  "\n"
  "Notice that there is no lazy evaluation nor optimization of the complex"
  " filters. Thus all parts of a complex filter are evaluated. For example,"
  " let us consider the filter expression:\n"
  "  '<expr1> and (<expr2> or <expr3>)'\n"
  "\n"
  "Both '<expr3>' and '<expr4>' are evaluated even if '<expr1>' is false. This"
  " may lead to waste of time and resources. Since multiple filters can be"
  " provided, it is more efficient to provide two simpler filters on command"
  " line:\n"
  "  --variant-filter '<expr1>' --variant-filter '<expr2> or <expr3>'\n"
  "\n"
  "Since all --variant-filter options must be satisfied, the variant filtering"
  " process is interrupted as soon as some filter discard a variant.\n"
  "\n"
  "In the end, even if the kim program makes possible to filter variants on the"
  " fly while creating the index, it should be more efficient to create the"
  " filtered VCF file(s) first by using any third party software of your choice,"
  " then use this(these) filtered VCF file(s) to create the index.\n"
  "\n"
  "Custom Variant Allele Frequencies\n"
  "---------------------------------\n"
  "\n"
  "You may want to compute identification metrics for various populations "
  "(e.g., \"EUR_AF\", \"AFR_AF\", ...) or use a VCF file that uses a non "
  "standard info tag to define the variant allele frequencies (e.g., "
  "\"FREQ\", ...). In such case, you can provide the VCF info tag to use "
  "in addition to the standard (\"AF\" [or \"AC\"/\"AN\"]) allele "
  "frequencies. The following example uses \"FR_FREQ\" and \"AU_FREQ\" to "
  "detect the French and Australian frequencies:\n"
  "\n"
  "\n  kim --create-index \\"
  "\n      --kmer-length 27 --kmer-prefix-length 6 \\"
  "\n      --reference file1.fasta --reference file2.fasta \\"
  "\n      --variants variants_file.vcf \\"
  "\n      --allele-frequency-tag FR_FREQ \\"
  "\n      --allele-frequency-tag AU_FREQ \\"
  "\n      --index-dir /where/to/store/index/"
};

const option::Descriptor KimProgram::_OPTION_FOOTER_QUERY = {
  UNKNOWN_OPT, 0, "" , "", Arg::None,
  "INDEX QUERY\n"
  "===========\n"
  "\n"
  "To query an existing index, the kim program needs an index and some"
  " files to analyse. An example of command line to query some index is:\n"
  "\n  kim --index-dir /path/to/my/index/ file1.fastq file2.fastq\n"
  "\n"
  "As for the index creation, the biological sequence files must be either"
  " fasta or fastq formatted.\n"
  "\n"
  "The default type I error (risk to reject a variant that really is in the data)"
  " is " KIM_DEFAULT_ALPHA " but it can be changed by using the --alpha (or -a)"
  " option:\n"
  "\n  kim --index-dir /path/to/my/index/ --alpha 0.001 file1.fastq file2.fastq\n"
  "\n  kim --index-dir /path/to/my/index/ --alpha 1e-3 file1.fastq file2.fastq\n"
  "\n  kim --index-dir /path/to/my/index/ --alpha 0.1% file1.fastq file2.fastq\n"
  "\n"
  "The default k-mer rate minimal threshold to consider a variant being in some"
  " read is " KIM_DEFAULT_THRESHOLD ". It can be changed too by using the"
  " --threshold (or -t) option:\n"
  "\n  kim --index-dir /path/to/my/index/ --threshold 0.3 file1.fastq file2.fastq\n"
  "\n  kim --index-dir /path/to/my/index/ --threshold 3e-1 file1.fastq file2.fastq\n"
  "\n  kim --index-dir /path/to/my/index/ --threshold 30% file1.fastq file2.fastq\n"
  "\n"
  "Using the strict mode can reduce the number of false positive but in "
  "the same time can increase the number of false negative. The for very "
  "small value of k (e.g., less than 20 for human genome), running in "
  "strict mode may improve the query results. But for \"standard\" values "
  "of k, running in weak mode (default) should give better results. The "
  "mode can be changed using:\n"
  "\n  kim --index-dir /path/to/my/index/ --mode strict file1.fastq file2.fastq\n"
};

const option::Descriptor KimProgram::_OPTION_END = {0, 0, 0, 0, 0, 0};

const option::Descriptor KimProgram::_short_usage[] = {
  _OPTION_PREAMBLE,
  _OPTION_INFORMATIONS_HEADER,
  _OPTION_INFORMATIONS_HELP,
  _OPTION_INFORMATIONS_VERSION,
  _OPTION_INFORMATIONS_COPYRIGHT,
  _OPTION_EMPTY_LINE,
  _OPTION_COMMON_OPTIONS_HEADER,
  _OPTION_COMMON_OPTIONS_QUIET,
  _OPTION_COMMON_OPTIONS_FORCE,
  _OPTION_COMMON_OPTIONS_INDEX_DIRECTORY,
  _OPTION_COMMON_OPTIONS_CHECK_CONSISTENCY,
  _OPTION_EMPTY_LINE,
  _OPTION_FOOTER_ABOUT,
  _OPTION_EMPTY_LINE,
  _OPTION_END
};

const option::Descriptor KimProgram::_create_index_usage[] = {
  _OPTION_PREAMBLE,
  _OPTION_INFORMATIONS_HEADER,
  _OPTION_INFORMATIONS_HELP,
  _OPTION_INFORMATIONS_VERSION,
  _OPTION_INFORMATIONS_COPYRIGHT,
  _OPTION_EMPTY_LINE,
  _OPTION_COMMON_OPTIONS_HEADER,
  _OPTION_COMMON_OPTIONS_QUIET,
  _OPTION_COMMON_OPTIONS_FORCE,
  _OPTION_COMMON_OPTIONS_INDEX_DIRECTORY,
  _OPTION_COMMON_OPTIONS_CHECK_CONSISTENCY,
  _OPTION_EMPTY_LINE,
  _OPTION_INDEX_CREATION_OPTIONS_HEADER,
  _OPTION_INDEX_CREATION_OPTIONS_CREATE_INDEX,
  _OPTION_INDEX_CREATION_OPTIONS_KMER_LENGTH,
  _OPTION_INDEX_CREATION_OPTIONS_KMER_PREFIX_LENGTH,
  _OPTION_INDEX_CREATION_OPTIONS_KMER_FILTER,
  _OPTION_INDEX_CREATION_OPTIONS_REFERENCE,
  _OPTION_INDEX_CREATION_OPTIONS_VARIANTS,
  _OPTION_INDEX_CREATION_OPTIONS_VARIANT_FILTER,
  _OPTION_INDEX_CREATION_OPTIONS_ALLELE_FREQUENCY_TAG,
  _OPTION_EMPTY_LINE,
  _OPTION_FOOTER_INDEX_CREATION,
  _OPTION_EMPTY_LINE,
  _OPTION_END
};

const option::Descriptor KimProgram::_query_usage[] = {
  _OPTION_PREAMBLE,
  _OPTION_INFORMATIONS_HEADER,
  _OPTION_INFORMATIONS_HELP,
  _OPTION_INFORMATIONS_VERSION,
  _OPTION_INFORMATIONS_COPYRIGHT,
  _OPTION_EMPTY_LINE,
  _OPTION_COMMON_OPTIONS_HEADER,
  _OPTION_COMMON_OPTIONS_QUIET,
  _OPTION_COMMON_OPTIONS_FORCE,
  _OPTION_COMMON_OPTIONS_INDEX_DIRECTORY,
  _OPTION_COMMON_OPTIONS_CHECK_CONSISTENCY,
  _OPTION_EMPTY_LINE,
  _OPTION_QUERY_OPTIONS_HEADER,
  _OPTION_QUERY_OPTIONS_OUTPUT_DIR,
  _OPTION_QUERY_OPTIONS_ALPHA,
  _OPTION_QUERY_OPTIONS_THRESHOLD,
  _OPTION_QUERY_OPTIONS_MODE,
  _OPTION_QUERY_OPTIONS_DUMP_READ_ANALYSIS,
  _OPTION_EMPTY_LINE,
  _OPTION_FOOTER_QUERY,
  _OPTION_EMPTY_LINE,
  _OPTION_END
};

const option::Descriptor KimProgram::_full_usage[] = {
  _OPTION_PREAMBLE,
  _OPTION_INFORMATIONS_HEADER,
  _OPTION_INFORMATIONS_HELP,
  _OPTION_INFORMATIONS_VERSION,
  _OPTION_INFORMATIONS_COPYRIGHT,
  _OPTION_EMPTY_LINE,
  _OPTION_COMMON_OPTIONS_HEADER,
  _OPTION_COMMON_OPTIONS_QUIET,
  _OPTION_COMMON_OPTIONS_FORCE,
  _OPTION_COMMON_OPTIONS_INDEX_DIRECTORY,
  _OPTION_COMMON_OPTIONS_CHECK_CONSISTENCY,
  _OPTION_EMPTY_LINE,
  _OPTION_INDEX_CREATION_OPTIONS_HEADER,
  _OPTION_INDEX_CREATION_OPTIONS_CREATE_INDEX,
  _OPTION_INDEX_CREATION_OPTIONS_KMER_LENGTH,
  _OPTION_INDEX_CREATION_OPTIONS_KMER_PREFIX_LENGTH,
  _OPTION_INDEX_CREATION_OPTIONS_KMER_FILTER,
  _OPTION_INDEX_CREATION_OPTIONS_REFERENCE,
  _OPTION_INDEX_CREATION_OPTIONS_VARIANTS,
  _OPTION_INDEX_CREATION_OPTIONS_VARIANT_FILTER,
  _OPTION_INDEX_CREATION_OPTIONS_ALLELE_FREQUENCY_TAG,
  _OPTION_EMPTY_LINE,
  _OPTION_QUERY_OPTIONS_HEADER,
  _OPTION_QUERY_OPTIONS_OUTPUT_DIR,
  _OPTION_QUERY_OPTIONS_ALPHA,
  _OPTION_QUERY_OPTIONS_THRESHOLD,
  _OPTION_QUERY_OPTIONS_MODE,
  _OPTION_QUERY_OPTIONS_DUMP_READ_ANALYSIS,
  _OPTION_EMPTY_LINE,
  _OPTION_FOOTER_ABOUT,
  _OPTION_EMPTY_LINE,
  _OPTION_FOOTER_INDEX_CREATION,
  _OPTION_EMPTY_LINE,
  _OPTION_FOOTER_QUERY,
  _OPTION_EMPTY_LINE,
  _OPTION_END
};

KimProgram::_OptionHandler::_OptionHandler(int argc, char **argv):
  stats(true, _full_usage, argc, argv),
  options(new option::Option[stats.options_max]),
  buffer(new option::Option[stats.buffer_max]),
  parse(_full_usage, argc, argv, options, buffer) {
}

KimProgram::_OptionHandler::~_OptionHandler() {
  delete [] options;
  delete [] buffer;
}


KimProgram::~KimProgram() {
}

KimProgram::KimProgram(int argc, char **argv):
  _progname(basename(argv[0])),
  _settings(), _running_mode(UNDEFINED_MODE),
  _kmer_filters(),
  _variant_filters(),
  _dna_files(), _variants_files(), _output_dir()
{
  _settings.setIndexDirectory(fs::weakly_canonical(KIM_DEFAULT_INDEX_DIRECTORY));

  argc -= (argc > 0);
  argv += (argc > 0); // skip program name argv[0] if present

  if (argc == 0) {
    option::printUsage(cout, _short_usage);
    exit(1);
  }

  _OptionHandler _opts(argc, argv);

  try {
    _processOptions(_opts);
  } catch (const exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }

}

string _lowercase(const char *s) {
  string res(s ? s : "");
  transform(res.begin(), res.end(), res.begin(), [](unsigned char c){ return tolower(c); });
  return res;
}

void KimProgram::_processOptions(_OptionHandler &_opts) {

  if (_opts.parse.error()) {
    throw Exception("Parse error");
  }

  if (_opts.options[HELP_OPT]) {
    string help_opt = _lowercase(_opts.options[HELP_OPT].arg);
    if (help_opt == "create-index") {
      option::printUsage(cout, _create_index_usage);
    } else if (help_opt == "query") {
      option::printUsage(cout, _query_usage);
    } else if (help_opt == "full") {
      option::printUsage(cout, _full_usage);
    } else {
      option::printUsage(cout, _short_usage);
    }
    exit(0);
  }

  if (_opts.options[VERSION_OPT]) {
    cout << "kim version " << VERSION << endl;
    exit(0);
  }

  if (_opts.options[COPYRIGHT_OPT]) {
    string copyright_opt = _lowercase(_opts.options[COPYRIGHT_OPT].arg);
    _showCopyright(copyright_opt == "full");
    exit(0);
  }

  if (_opts.options[UNKNOWN_OPT]) {
    Exception e("The kim program doesn't accept the following options:");
    for (option::Option* opt = _opts.options[UNKNOWN_OPT]; opt; opt = opt->next()) {
      e << " ";
      /*      if (opt->name[0] != '-') {
        e << "-";
        }*/
      e << opt->name;
    }
    throw e;
  }

  _settings.warn(!_opts.options[QUIET_OPT]);

  _settings.allowOverwrite(_opts.options[FORCE_OPT]);

  if (_settings.warn()) {
    cerr << _progname << " version " VERSION << endl
         << "Starting time: " << timestamp() << endl;
  }

  _settings.checkConsistency(_opts.options[CONSISTENCY_CHECKING_OPT]);
  if (_settings.checkConsistency() && _settings.warn()) {
    cerr << "Consistency Checking while reading DNA file is enabled. "
         << "This will significantly impact the program running time."
         << endl;
  }

  if (_opts.options[INDEX_DIRECTORY_OPT]) {
    _settings.setIndexDirectory(fs::weakly_canonical(_opts.options[INDEX_DIRECTORY_OPT].arg),
                                !_opts.options[CREATE_INDEX_OPT],
                                _opts.options[CREATE_INDEX_OPT] && !_settings.allowOverwrite());
  }

  if (_opts.options[CREATE_INDEX_OPT]) {

    _running_mode = INDEX_CREATION_MODE;

    // Ensure that at least one reference file is provided
    if (_opts.options[REFERENCE_OPT].count() < 1) {
      throw Exception("At least one reference file must be provided for index creation.");
    }

    // Ensure that at least one variant file is provided
    if (_opts.options[VARIANTS_OPT].count() < 1) {
      throw Exception("At least one variant file must be provided for index creation.");
    }

    // Add extra allele frequency tags
    _settings.addAlleleFrequencyTag("AF");
    for (option::Option* opt = _opts.options[ALLELE_FREQUENCY_TAG_OPT]; opt; opt = opt->next()) {
      if (!_settings.addAlleleFrequencyTag(opt->arg) && _settings.warn()) {
        cerr << "WARNING: Allele frequency tag '" << opt->arg << "' was probably given twice."
             << endl;
      }
    }

    // Reserve memory for k-mer exclusion filters, then fill the
    // _kmer_filters vector.
    _kmer_filters.reserve(_opts.options[KMER_FILTER_OPT].count());
    for (option::Option* opt = _opts.options[KMER_FILTER_OPT]; opt; opt = opt->next()) {
      regex filter(opt->arg, regex::optimize);
      _kmer_filters.emplace_back(pair(opt->arg, filter));
    }

    // Reserve memory for variant filters.
    _variant_filters.reserve(_opts.options[VARIANT_FILTER_OPT].count());

    // Check if all given variant fitler are valid then fill the
    // _variant_filters vector.
    vcfpp::BcfWriter bcf("/dev/null", "VCF4.3");
    bcf.header.addFORMAT("GT", "1", "String", "Genotype");
    bcf.header.addINFO("AF", "A", "Float", "Estimated allele frequency in the range (0,1)");
    bcf.header.addContig("A"); // add chromosome
    vcfpp::BcfRecord r(bcf.header);
    r.setCHR("A");
    r.setPOS(2);
    r.setID("fake");
    r.setRefAlt("A,C");
    r.setQUAL('@');
    r.setFILTER("PASS");
    for (option::Option* opt = _opts.options[VARIANT_FILTER_OPT]; opt; opt = opt->next()) {
      VariantFilterDriver drv(r);
      string filter = opt->arg;
      drv.apply(filter);// If filter is not syntaxically correct, an exception is thrown.
      _variant_filters.push_back(filter);
    }

    // Check if all given reference files are readable and fill the
    // dna_files vector with their filenames.
    _dna_files.reserve(_opts.options[REFERENCE_OPT].count());
    for (option::Option* opt = _opts.options[REFERENCE_OPT]; opt; opt = opt->next()) {
      fs::path p = opt->arg;
      _assertFileIsReadable(p);
      _dna_files.push_back(fs::canonical(p));
    }

    // Check if all given variant files are readable and fill the
    // variant_files vector with their filenames.
    _variants_files.reserve(_opts.options[VARIANTS_OPT].count());
    for (option::Option* opt = _opts.options[VARIANTS_OPT]; opt; opt = opt->next()) {
      if (!_isURL(opt->arg)) _assertFileIsReadable(opt->arg);
      _variants_files.push_back(opt->arg);
    }

    // Set the k-mer length
    size_t kmer_length = KIM_DEFAULT_KMER_LENGTH;
    if (_opts.options[KMER_LENGTH_OPT]) {
      kmer_length = stoull(_opts.options[KMER_LENGTH_OPT].arg);
    }
    _settings.setKmerLength(kmer_length);

    // Set the k-mer prefix length
    size_t kmer_prefix_length = KIM_DEFAULT_KMER_PREFIX_LENGTH;
    if (_opts.options[KMER_PREFIX_LENGTH_OPT]) {
      kmer_prefix_length = stoull(_opts.options[KMER_PREFIX_LENGTH_OPT].arg);
    }
    _settings.setKmerPrefixLength(kmer_prefix_length);

    // Freeze settings
    _settings.freeze();

    // Prevent using options defined for other mode(s) than index
    // creation.
    for (auto o: {
        OUTPUT_OPT,
        ALPHA_OPT,
        THRESHOLD_OPT,
        MODE_OPT,
        DUMP_READ_ANALYSIS_OPT,
      }) {
      if (_opts.options[o]) {
        Exception e;
        e << "Option '" << _opts.options[o].name << "'"
          " is not available for index creation.";
        throw e;
      }
    }

  } else {

    _running_mode = QUERY_MODE;

    // If not in index creation mode, then must be in query mode and
    // thus must have some file to process.
    if (_opts.parse.nonOptionsCount() == 0) {
      throw Exception("Bad usage");
    }

    // Check if all given query files are readable and fill the
    // dna_files vector with their filenames.
    _dna_files.reserve(_opts.parse.nonOptionsCount());
    for (int i = 0; i < _opts.parse.nonOptionsCount(); ++i) {
      fs::path p = _opts.parse.nonOption(i);
      _assertFileIsReadable(p);
      _dna_files.push_back(fs::canonical(p));
    }

    // Ensure that given output directory (if provided) is writable.
    if (_opts.options[OUTPUT_OPT]) {
      _output_dir = _opts.options[OUTPUT_OPT].arg;
      fs::create_directories(_output_dir);
      Settings::validateDirectory(_output_dir, true);
      _output_dir = fs::canonical(_output_dir);
    }

    string alpha_str = (_opts.options[ALPHA_OPT]
                        ? _opts.options[ALPHA_OPT].arg
                        : KIM_DEFAULT_ALPHA);
    double alpha = stod(alpha_str);
    if (alpha_str.back() == '%') {
      alpha *= 0.01;
    }
    _settings.alpha(alpha);

    string threshold_str = (_opts.options[THRESHOLD_OPT]
                            ? _opts.options[THRESHOLD_OPT].arg
                            : KIM_DEFAULT_THRESHOLD);
    double threshold = stod(threshold_str);
    if (threshold_str.back() == '%') {
      threshold *= 0.01;
    }
    _settings.threshold(threshold);

    string mode_str = (_opts.options[MODE_OPT]
                       ? _lowercase(_opts.options[MODE_OPT].arg)
                       : KIM_DEFAULT_MODE);
    if (mode_str == "weak") {
      _settings.weakMode(true);
    } else if (mode_str == "strict") {
      _settings.strictMode(true);
    } else {
      Exception e;
      e << "Incorrect value '" << mode_str << "' "
        << "for option '" << _opts.options[MODE_OPT].name << "'.";
      throw e;
    }

    _dump_read_analysis_opt = _opts.options[DUMP_READ_ANALYSIS_OPT];

    // Prevent using options defined for other mode(s) than index
    // querying.
    for (auto o: {
        VARIANT_FILTER_OPT,
        KMER_LENGTH_OPT,
        KMER_PREFIX_LENGTH_OPT,
        KMER_FILTER_OPT,
        REFERENCE_OPT,
        VARIANTS_OPT,
        ALLELE_FREQUENCY_TAG_OPT
      }) {
      if (_opts.options[o]) {
        Exception e;
        e << "Option '" << _opts.options[o].name << "'"
          " is not available for index querying.";
        throw e;
      }
    }

  }

}

bool KimProgram::_variantOfInterest(vcfpp::BcfRecord &v) const {
  bool res = true;
  vector<string>::const_iterator it = _variant_filters.begin();
  VariantFilterDriver drv(v);
  while (res && it != _variant_filters.end()) {
    res = drv.apply(*it);
    ++it;
  }
  return res;
}

bool KimProgram::_isURL(const string &f) {
  static const regex url_regex("^(([^:/?#]+):)(//([^/?#]*))?([^?#]*)(\\?([^#]*))?(#(.*))?$"); // From https://www.rfc-editor.org/rfc/rfc2396#appendix-B
  return regex_match(f, url_regex);
}

bool KimProgram::_checkIfFileIsReadable(const fs::path &f) {
  bool res = false;
  ifstream ifs(f);
  if (ifs) {
    ifs.close();
    res = true;
  }
  return res;
}

void KimProgram::_assertFileIsReadable(const fs::path &f) {
  if (!_checkIfFileIsReadable(f)) {
    Exception e;
    e << "Unable to find or open the file " << f << ".";
    throw e;
  }
}

// Dump the file content (up to the number of given characters) to the
// given stream.
void KimProgram::_dumpFile(const fs::path &filename, size_t n, ostream &os) {
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

void KimProgram::_dumpResultsHeader(const string &input_fname, ostream &os) const {
  os << "---\n"
     << "Informations:\n"
     << "  File: " << input_fname << "\n"
     << "  Index: " << _settings.getIndexDirectory() << "\n"
     << "  k-mer length: " << _settings.k() << "\n"
     << "  Alpha: " << _settings.alpha() << "\n"
     << "  Threshold: " << _settings.threshold() << "\n"
     << "  Mode: " << (_settings.weakMode() ? "Weak" : "Strict") << "\n";
}

void KimProgram::_dumpResultsRead(const ReadAnalyzer::VariantKmerRates &res, const _currentReadInfos &infos, ostream &os) const {
  os << "- Read: " << infos.description << "\n"
     << "  Id: " << infos.id << "\n"
     << "  Variants:\n";
  for (ReadAnalyzer::VariantKmerRates::const_iterator it = res.cbegin();
       it != res.end();
       ++it) {
    ReadAnalyzer::VariantKmerRatesConstIteratorWrapper vs(it);
    os << "  - id: "  << vs.node.variant << "\n"
       << "    associated k-mers: " << vs.node.in_degree << "\n"
       << "    found k-mer rate: " << vs.rate << "\n"
       << "    z-score: " << vs.zscore << "\n";
  }
}

double log2(double x) {
  return log(x) / log(2);
}

void KimProgram::_dumpResultsVariants(const ReadAnalyzer::VariantScores &results, ostream &os) const {

  os << "Variants:" << (results.empty() ? " []" : "") << "\n";

  /*
   * Let consider the Bernouilli distribution of having $X$ k-mers
   * associated to the current variant in some read due to random.
   *
   * The parameter of this distribution are:
   * - $n$ The number of k-mers in a read
   * - $p$ The probability of some k-mer to be one of the variant
   *   associated k-mer
   * - $X$ The variable following a binomial distribution of
   *   parameters n and p (the number of indexed k-mers in a read
   *   having $n$ k-mers due to random).
   *
   * Thus $E[X] = n\,p$ and $V[X] = n\,p\,(1-p)$.
   *
   * The Bienaymé–Chebyshev inequality states that:
   * \begin{align*}
   *   P\left(|X - n\,p| > x\right)                  & \leq \frac{n\,p\,(1 - p)}{{(n\,x)}^2} \\
   *   P\left(|X - n\,p| > x\right)                  & \leq \frac{n\,p\,(1 - p)}{{(n\,x)}^2} \\
   *   P\left(|\frac{X}{n} - p| > \frac{x}{n}\right) & \leq \frac{n\,p\,(1 - p)}{{(n\,x)}^2} \\
   *   P\left(|\frac{X}{n} - p| > \frac{x}{n}\right) & \leq \frac{p\,(1 - p)}{n\,x^2} \\
   * \end{align*}
   *
   * Here $x$ is the observed number of associated k-mers, thus
   * \frac{x}{n} is the observed rate of associated k-mers due to
   * random in the read.
   *
   * For convenience, let define $Y = \frac{X}{n}$ and $y
   * =\frac{x}{n}$ (and recall that $n \geq k$). This gives:
   * \begin{align*}
   *               P\left(|Y - p| > y\right) & \leq \frac{p\,(1 - p)}{n^3\,y^2} \\
   *   \Rightarrow P\left(|Y - p| > y\right) & \ll \frac{p\,(1 - p)}{k^3\,y^2}
   * \end{align*}
   *
   * The $\frac{p\,(1 - p)}{k^3\,y^2}$ value is an over-estimate
   * upper bound of the probability of seeing the k-mer associated
   * to the current variant due to random, thus the under-estimate
   * confidence lower bound is its complement to 1.
   */

  static const size_t nb_kmers = 1 << (_settings.k() << 1);

  for (ReadAnalyzer::VariantScores::const_iterator it = results.cbegin();
       it != results.cend();
       ++it) {
    ReadAnalyzer::VariantScoreConstIteratorWrapper vs(it);
    double p = double(vs.node.in_degree) / nb_kmers;
    double uncertainty = (p * (1.0 - p) / vs.stats.mean / vs.stats.mean / _settings.k() / _settings.k() / _settings.k());
    os << "- id: " << vs.node.variant << "\n"
       << "    associated k-mers: " << vs.node.in_degree << "\n"
       << "  stats:" << "\n"
       << "    mean: " << vs.stats.mean << "\n"
       << "    variance: " << vs.stats.variance << "\n"
       << "    count: " << vs.stats.count << "\n"
       << "    uncertainty: " << uncertainty << "\n"
      ;
  }
}

void KimProgram::_dumpResultsKimScores(const map<string, double> &scores, ostream &os) {
  os << "Identification metrics:\n";
  for (auto const &[tag, score]: scores) {
    if (!isnan(score)) {
      os << "- " << tag << ":\n"
         << "    probability: " << score << "\n"
         << "    e-value (1G. people): " << (score * 1e9) << "\n";
    }
  }
}

void KimProgram::_dumpResultsFooter(Monitor &monitor, ostream &os) {
  os << "Monitoring:" << "\n"
     << "  Wallclock time: "
     << (monitor.getWallClockTime() / 1s) << "s" << "\n"
     << "  User CPU time: " << (monitor.getUserTime() / 1s) << "s" << "\n"
     << "  System CPU time: " << (monitor.getSystemTime() / 1s) << "s" << endl;
}

void KimProgram::_showCopyright(bool full) {
   cout << PACKAGE " version " VERSION " -- " PACKAGE_SHORT_DESCRIPTION "\n\n";
   string fname_orig = "LICENSE.md";
   string fname = FileReader::findFile(fname_orig);
  if (!fname.empty()) {
    _dumpFile(fname, 508);
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
      _dumpFile(fname);
    } else {
      cerr << "File '" << fname_orig << "' not found or not readable."
           << " Please ensure the kim program is correctly installed."
           << endl;
    }
    cout << endl;
    fname_orig = "LICENSE-fr.md";
    fname = FileReader::findFile(fname_orig);
    if (!fname.empty()) {
      _dumpFile(fname);
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

fs::path KimProgram::_commonPathPrefix(const vector<string> &files) {
  fs::path p;
  if (!files.empty()) {
    p = fs::path(files[0]).parent_path();
    for (size_t i = 1; i < files.size(); ++i) {
      fs::path q = fs::path(files[i]).parent_path();
      fs::path::iterator p_it = p.begin();
      fs::path::iterator q_it = q.begin();
      fs::path new_p;
      while ((p_it != p.end()) && (q_it != q.end()) && (*p_it == *q_it)) {
        new_p /= *p_it;
        ++p_it;
        ++q_it;
      }
      p = new_p;
    }
  }
  return p;
}

fs::path KimProgram::_buildResultFile(const fs::path &f) const {
  assert(!f.empty());
  assert(f.is_absolute());
  fs::path p = f.relative_path();
  p.replace_extension(KIM_DEFAULT_RESULT_EXTENSION);
  fs::path out_file = _output_dir / "";
  string sep = "";
  for (fs::path::iterator it = p.begin(); it != p.end(); ++it) {
    out_file += sep;
    out_file += *it;
    sep = "_";
  }

  if (_checkIfFileIsReadable(out_file)) {
    if (_settings.allowOverwrite()) {
#ifdef _OPENMP
#  pragma omp critical
#endif
      cerr << "File " << out_file << " already exist and will be overwritten." << endl;
    } else {
      Exception e;
      e << "File " << out_file << " already exist. You should either provide a new output directory using the --output-dir option (recommanded) or use the --force flag (dangerous)" << endl;
      throw e;
    }
  }
  return out_file;
}

void _updateReadInfosCallback(const DNAFileReader &reader, _currentReadInfos &infos) {
  infos.description = reader.getCurrentSequenceDescription();
  infos.id = reader.getCurrentSequenceID();
}

void _readAnalyzerCallback(ReadAnalyzer &r) {
  r.analyze();
  r.reset();
}

void _readAnalyzerDumpCallback(const KimProgram &prog, ReadAnalyzer &r, _currentReadInfos &infos, ostream &os) {
  const ReadAnalyzer::VariantKmerRates &res = r.analyze();
  if (!res.empty()) {
    prog._dumpResultsRead(res, infos, os);
    r.reset();
  }
}

void KimProgram::_runQuery() {
  Monitor monitor;
  // Loading the index.
  KmerVariantGraph kim_index(_settings);
  size_t common_path_prefix_size = _commonPathPrefix(_dna_files).string().size();

  kim_index.freeze();

  double indexed_kmer_probability = kim_index.getNbKmers();
  indexed_kmer_probability /= (1 << (_settings.k() << 1));
  // Process each input file
#ifdef _OPENMP
  if (_output_dir.empty()) {
    // Disable multithreading if results are output on the terminal
    omp_set_num_threads(1);
  }
#  pragma omp parallel for
#endif
  for (auto const &f: _dna_files) {

    ofstream ofs;
    // Write the result to output stream
    if (!_output_dir.empty()) {
      fs::path out_file = _buildResultFile(f.substr(common_path_prefix_size));
      ofs.open(out_file);
      if (!ofs) {
        Exception e;
        e << "Unable to open " << _output_dir << " file to print results.";
        throw e;
      }
    }
    ostream &os = ofs.is_open() ? ofs : cout;

    _dumpResultsHeader(f, os);

    Monitor one_file_monitor;

    if (_settings.warn()) {
#ifdef _OPENMP
#  pragma omp critical
#endif
      cerr << "Processing file '" << f << "'" << endl;
    }
    DNAFileReader reader(_settings.k(), f, _settings.warn());
    ReadAnalyzer read_analyzer(_settings.alpha(), _settings.threshold(), _settings.weakMode());
    _currentReadInfos read_infos;

    reader.gotoNextSequence();
    read_infos.description = reader.getCurrentSequenceDescription();
    read_infos.id = reader.getCurrentSequenceID();

    const DNAFileReader::onNewSequenceFct readAnalyzerFct = [&read_analyzer](const DNAFileReader &) {
      _readAnalyzerCallback(read_analyzer);
    };

    const DNAFileReader::onNewSequenceFct readAnalyzerDumpFct = [this, &read_analyzer, &read_infos, &os](const DNAFileReader &) {
      _readAnalyzerDumpCallback(*this, read_analyzer, read_infos, os);
    };
    const DNAFileReader::onNewSequenceFct updateReadInfosFct = [&read_infos](const DNAFileReader &reader) {
      _updateReadInfosCallback(reader, read_infos);
    };

    size_t pass_number = 0, nb_pass = _dump_read_analysis_opt ? 2 : 1;

    do {
      reader.reset();
      reader.removeOnNewSequenceCallback(readAnalyzerFct);
      reader.removeOnNewSequenceCallback(readAnalyzerDumpFct);
      reader.removeOnNewSequenceCallback(updateReadInfosFct);
      switch (++pass_number) {
      case 1:
        reader.addOnNewSequenceCallback(readAnalyzerFct);
        break;
      case 2:
        os << "Read Analysis:\n";
        reader.addOnNewSequenceCallback(readAnalyzerDumpFct);
        reader.addOnNewSequenceCallback(updateReadInfosFct);
        break;
      default:
        throw Exception("BUG");
      }

      // Process each k-mer from the current input file
      for (string kmer = reader.getNextKmer(true /* skip degenerate */);
           !kmer.empty();
           kmer = reader.getNextKmer(true /* skip degenerate */)) {
        // cerr << "k-mer : '" << kmer << "'" << endl;
        list<KmerVariantEdgesSubindex::KmerVariantAssociation> assoc_variant = kim_index.search(kmer);
        // cerr << "assoc_variant size is " << assoc_variant.size() << endl;

        // We don't care of in_ref if the k-mer is not associated to som
        // variant, thus we can assign any value. Although, when the
        // k-mer is not indexed, calling isInReferenceKmer() will throw
        // an exception. So the lazy evaluation when assigning true will
        // prevent this exception.
        bool in_ref = assoc_variant.empty() || kim_index.isInReferenceKmer(kmer);
        for (list<KmerVariantEdgesSubindex::KmerVariantAssociation>::const_iterator it = assoc_variant.cbegin(); it != assoc_variant.cend(); ++it) {
          // cerr << "Current variant association concerns variant '" << it->variant_node.variant << "' which concerns " << it->variant_node.in_degree << " k-mers" << endl;
          // cerr << "This variant is potentially seen in read " << reader.getCurrentSequenceDescription()
          //      << " by the k-mer at position " << reader.getCurrentKmerRelativePosition()
          //      << " in the read" << endl;
          read_analyzer.add(it->variant_node, it->rank, in_ref);
        }
      }

      switch (pass_number) {
      case 1:
        _readAnalyzerCallback(read_analyzer);
        read_analyzer.complete();
        break;
      case 2:
        _readAnalyzerDumpCallback(*this, read_analyzer, read_infos, os);
        _updateReadInfosCallback(reader, read_infos);
        break;
      default:
        throw Exception("BUG");
      }
    } while (pass_number < nb_pass);

    one_file_monitor.stop();

    // Analyse reads for variants
    const ReadAnalyzer::VariantScores &result = read_analyzer.result();

    // Compute the kim scores
    map<string, double> scores = _computeKimScores(kim_index, result);

    // Write the result and scores to output stream
    _dumpResultsVariants(result, os);
    _dumpResultsKimScores(scores, os);
    _dumpResultsFooter(one_file_monitor, os);

    if (ofs.is_open()) {
      ofs.close();
    }
  }

  monitor.stop();
  cerr << "# Resources used for querying the index:" << endl
       << "# - Wallclock time (in seconds): "
       << (monitor.getWallClockTime() / 1s) << endl
       << "# - User CPU time (in seconds): " << (monitor.getUserTime() / 1s) << endl
       << "# - System CPU time (in seconds): " << (monitor.getSystemTime() / 1s) << endl
       << "# - Memory: " << Monitor::memoryWithUnit2string(monitor.getMemory()) << endl;
}

void KimProgram::_createTemporaryDirectory() {
  const string tmpdir_template = (fs::temp_directory_path() / "kimAF_XXXXXX").string();
  char *tmpdir = new char[tmpdir_template.size() + 1];
  copy(tmpdir_template.begin(), tmpdir_template.end(), tmpdir);
  tmpdir[tmpdir_template.size()] = '\0';
  if (mkdtemp(tmpdir) == NULL) {
    Exception e;
    e << "The following error occurs while trying to create the temporary directory '" << tmpdir << "': " << strerror(errno);
    delete [] tmpdir;
    throw e;
  }
  _tmpdir = tmpdir;
  delete [] tmpdir;
  if (_settings.warn()) {
    cerr << "* Temporary directory " << _tmpdir << " created." << endl;
  }
}

void KimProgram::_deleteTemporaryDirectory() {
  if (!fs::remove(_tmpdir)) {
    Exception e;
    e << "Unable to remove the directory '" << _tmpdir << ".";
    throw e;
  }
  if (_settings.warn()) {
    cerr << "* Temporary directory " << _tmpdir << " deleted." << endl;
  }
}

void KimProgram::_createTemporaryVafFiles(KmerVariantGraph &index) {
  bool frozen = index.frozen();
  if (frozen) index.unfreeze();
  if (_tmpdir.empty()) {
    _createTemporaryDirectory();
  }
  assert(!_tmpdir.empty());
  fs::path current_file = _tmpdir / VARIANT_ALLELE_FREQUENCIES_BASENAME;
  assert(_vaf_files.find("") == _vaf_files.end());
  fstream &os = _vaf_files[""];
  os.open(current_file, ios::out | ios::trunc);
  if (!os.is_open()) {
    Exception e;
    e << "An error occurs while trying to create the file " << current_file << ".";
    throw e;
  }
  index.addVariantAlleleFrequencyFile("", VARIANT_ALLELE_FREQUENCIES_BASENAME);

  for (auto const &tag: _settings.getAlleleFrequencyTags()) {
    fs::path current_file_basename = _buildVafFilename(tag);
    current_file = _tmpdir / current_file_basename;
    assert(_vaf_files.find(tag) == _vaf_files.end());
    fstream &os = _vaf_files[tag];
    os.open(current_file, ios::out | ios::trunc | ios::binary);
    if (!os.is_open()) {
      _closeVafFiles();
      Exception e;
      e << "An error occurs while trying to create the file " << current_file << ".";
      throw e;
    }
    index.addVariantAlleleFrequencyFile(tag, current_file_basename);
  }
  if (frozen) index.freeze();
}

void moveFile(const fs::path &old_p, const fs::path &new_p) {
  try {
    fs::rename(old_p, new_p);
  } catch (...) {
    fs::copy_file(old_p, new_p);
    fs::remove(old_p);
  }
}

void KimProgram::_saveTemporaryVafFiles() {
  assert(!_tmpdir.empty());
  assert(!_settings.getIndexDirectory().empty());
  if (_settings.warn()) {
    cerr << "* Saving variant allele frequencies to '" << _settings.getIndexDirectory() << "'." << endl;
  }
  _closeVafFiles();
  fs::path current_file = VARIANT_ALLELE_FREQUENCIES_BASENAME;
  moveFile(_tmpdir / current_file, _settings.getIndexDirectory() / current_file);
  for (auto const &tag: _settings.getAlleleFrequencyTags()) {
    current_file = _buildVafFilename(tag);
    moveFile(_tmpdir / current_file, _settings.getIndexDirectory() / current_file);
    if (_settings.warn()) {
      cerr << "  - File " << current_file << " saved." << endl;
    }
  }
}

void KimProgram::_loadIndexVafFiles(const KmerVariantGraph &index) {
  assert(!_settings.getIndexDirectory().empty());
  fs::path _index_dir = _settings.getIndexDirectory();
  for (auto const &[tag, file]: index.variantAlleleFrequencyFiles()) {
#ifdef DEBUG
    cerr << "Opening file " << file << " for " << (tag.empty() ? "Variant List" : tag) << endl;
#endif
    string t = tag.empty() ? VARIANT_ALLELE_FREQUENCIES_BASENAME : tag;
    fs::path current_file = _index_dir / file;
    fstream &is = _vaf_files[tag];
    assert(!is.is_open());
    is.open(current_file, ios::in | ios::binary);
    if (!is.is_open()) {
      Exception e;
      e << "An error occurs while trying to open the file " << current_file << ".";
      throw e;
    }
  }
}

void KimProgram::_closeVafFiles() {
  for (auto &[fname, fstr]: _vaf_files) {
#ifdef DEBUG
    cerr << "Closing file " << (fname.empty() ? "Variant List": fname) << endl;
#endif
    fstr.close();
  }
}

static const uint16_t encoded_af_delta = 32768;
static const uint16_t encoded_af_missing = 32767;
static const uint16_t encoded_af_max = 32766;

void KimProgram::_dump2temporaryVafFiles(VariantAlleleFrequencies &vaf) {
  if (vaf.getCurrentVariantAlleleFrequencies().empty()) return;
  assert(_vaf_files.find("") != _vaf_files.end());
  fstream &v_os = _vaf_files[""];
  assert(v_os.is_open());
  uint16_t delta = 0;
  for (auto const &vpaf: vaf.getCurrentVariantAlleleFrequencies()) {
    if (vpaf.variant != "-") {
      v_os << vpaf.variant << "\n";
      for (auto const &tag: _settings.getAlleleFrequencyTags()) {
        assert(_vaf_files.find(tag) != _vaf_files.end());
        fstream &os = _vaf_files[tag];
        assert(os.is_open());
        map<string, float>::const_iterator it = vpaf.population_af.find(tag);
        uint16_t encoded_af = ((it == vpaf.population_af.cend()) ? encoded_af_missing : (it->second * encoded_af_max));
        encoded_af += delta;
        unsigned char bytes[2];
	bytes[0] = (unsigned char)(encoded_af >> 8);
	bytes[1] = (unsigned char)(encoded_af);
        os.write(reinterpret_cast<char*>(bytes), 2);
      }
      delta = encoded_af_delta;
    }
  }
}

double __extractFrequency(istream &is) {
  assert(is);
  unsigned char bytes[2];
#ifdef DEBUG
  cerr << "Reading two bytes" << endl;
#endif
  is.read(reinterpret_cast<char*>(bytes), 2);
  assert(is);
  uint16_t encoded_af = (uint16_t)((bytes[0] << 8) | bytes[1]);
#ifdef DEBUG
  cerr << "Value is " << (int) encoded_af << endl;
#endif
  if (encoded_af >= encoded_af_delta) encoded_af -= encoded_af_delta;
  if (encoded_af == encoded_af_missing) return -1;
  double af = double(encoded_af) / double(encoded_af_max);
#ifdef DEBUG
  cerr << "Decoded value is " << af << endl;
#endif
  assert(af >= 0);
  assert(af <= 1.);
  return af;
}

map<string, double> KimProgram::_computeKimScores(const KmerVariantGraph &index, const ReadAnalyzer::VariantScores &result) {
  map<string, double> scores;
  for (auto const &tag: _settings.getAlleleFrequencyTags()) {
#ifdef DEBUG
    cerr << "Initialization of scores[" << tag << "] to 1" << endl;
#endif
    scores[tag] = nan("");
  }
  _loadIndexVafFiles(index);
  assert(_vaf_files.find("") != _vaf_files.cend());
  istream &v_is = _vaf_files[""];
#ifdef DEBUG
  cerr << "Parsing " << VARIANT_ALLELE_FREQUENCIES_BASENAME << " file" << endl;
#endif
  while (v_is) {

    // Using "AF" tag to detect alternative alleles.
    assert(_vaf_files.find("AF") != _vaf_files.cend());
    fstream &af_is = _vaf_files["AF"];

    size_t nb_alt = -1;
    size_t variant_pos = -1;
    string variant;

    double af = 0.0;
    double af_complement = 1.0;

    do {
      ++nb_alt;
      string line;
      getline(v_is, line);
      if (v_is) {
#ifdef DEBUG
        cerr << "Checking variant '" << line << "'" << endl;
#endif
        assert(v_is);
        assert(!line.empty());
        double cur_af = __extractFrequency(af_is);
        if (cur_af >= 0) {
          assert(cur_af <= af_complement);
          af_complement -= cur_af;
#ifdef DEBUG
          cerr << "AF frequency is " << cur_af << " and complement is now " << af_complement << endl;
#endif
          const VariantNodesIndex::const_iterator &it = index.getVariantNodesIndex().find(line);
          if ((it != index.getVariantNodesIndex().cend()) && (result.find(it) != result.end())) {
#ifdef DEBUG
            cerr << "*** Variant '" << line << "' is at rank " << nb_alt << endl;
#endif
            variant_pos = nb_alt;
            assert(variant.empty());
            variant = line;
            assert(!variant.empty());
            af = cur_af;
#ifdef DEBUG
            cerr << "The AF is " << af << endl;
#endif
          }
        }
      }
    } while (v_is && af_is && (af_is.peek() & 128));

    assert(!v_is || (nb_alt != (size_t) -1));
    if (v_is) {
#ifdef DEBUG
      cerr << "The score for AF is " << scores["AF"] << endl;
      cerr << "variant is " << variant << " (at pos " << variant_pos << ")" << endl;
      cerr << "af_complement is " << af_complement << " and af " << af << endl;
#endif
      if (af >= 0) {
        double &s = scores["AF"];
        if (isnan(s)) s = 1.0;
        s *= (variant.empty() ? af_complement : af);
      }

      // Compute the probability for the other tags
      for (auto const &tag: _settings.getAlleleFrequencyTags()) {
        if (tag != "AF") {
          af = -1.0;
          af_complement = 1.0;
          for (size_t i = 0; i < nb_alt; ++i) {
#ifdef DEBUG
            cerr << "Checking af for '" << tag << "'" << endl;
#endif
            assert(_vaf_files.find(tag) != _vaf_files.cend());
            fstream &tag_is = _vaf_files[tag];
            double cur_af = __extractFrequency(tag_is);
            if (cur_af >= 0) {
              assert(cur_af <= af_complement);
              af_complement -= cur_af;
#ifdef DEBUG
              cerr << tag << " frequency is " << cur_af << " and complement is now " << af_complement << endl;
#endif
              if (i == variant_pos) {
                af = cur_af;
#ifdef DEBUG
                cerr << "The " << tag << " is " << af << endl;
#endif
              }
            }
          }
          if (af >= 0) {
            double &s = scores[tag];
            if (isnan(s)) s = 1.0;
            s *= (variant.empty() ? af_complement : af);
          }
        }
#ifdef DEBUG
        cerr << "The score for " << tag << " is now " << scores[tag] << endl;
#endif
      }
    }
  }
  _closeVafFiles();
  return scores;
}


void KimProgram::_createIndex() {

  Monitor monitor;

  // Checks validity of current settings
  Settings fake_settings(_settings);
  try {
    fake_settings.unfreeze();
    fake_settings.setIndexDirectory(_settings.getIndexDirectory(), false, true);
  } catch (BadSettingsException &e) {
    // The index directory already exists.
    if (_settings.allowOverwrite()) {
      fake_settings.setIndexDirectory(_settings.getIndexDirectory()
                                      / "this kind of path should not exist",
                                      false, true);
      KmerVariantGraph fake_index(fake_settings);
      if (_settings.warn()) {
        cerr << "* ";
      }
      fake_index.removeDumpedIndex(_settings.getIndexDirectory());
      fake_settings.unfreeze();
      fake_settings.setIndexDirectory(_settings.getIndexDirectory(), false, true);
    } else {
      throw;
    }
  }

  // Prepare variant k-mer enumaration
  VariantKmerEnumerator::init(_settings, _dna_files);

  // Create an empty k-mer to variant graph
  KmerVariantGraph kim_index(_settings);
  kim_index.freeze();
  assert(kim_index.frozen());
  assert(kim_index.getNbKmerVariantEdges() == 0);
  assert(kim_index.getNbKmers() == 0);
  assert(kim_index.getNbVariants() == 0);

  // Prepare variant allele frequency computation
  VariantAlleleFrequencies vaf(_settings.getAlleleFrequencyTags());
  _createTemporaryVafFiles(kim_index);

  // Processing Variants files
  for (auto const &f: _variants_files) {
    if (_settings.warn()) {
      cerr << "* Processing variant file '" << f << "'" << endl;
    }
    vcfpp::BcfReader vcf(f);
#ifdef _OPENMP
    vcf.setThreads(omp_get_num_threads());
#endif
    const vcfpp::BcfHeader &hdr = vcf.header;

    // set<string> usable_tags;
    int tag_type = hdr.getInfoType("AF");
    bool allele_frequency_available = (tag_type != -1)
      || ((hdr.getInfoType("AC") != -1) && (hdr.getInfoType("AN") != 1));
    bool allele_frequency_usable = (tag_type == 2)
      || ((hdr.getInfoType("AC") == 1) && (hdr.getInfoType("AN") == 1));
    // if (allele_frequency_usable) {
    //   usable_tags.insert("AF");
    // }
    if (_settings.warn()) {
      cerr << "  - Variant Allele frequency:\n";
      if (allele_frequency_available && !allele_frequency_usable) {
        cerr << "WARNING: The default allele frequency doesn't follows the VCF specification and thus can't be used.\n";
      }
      cerr << "  - Variant Allele frequency:\n"
           << "    - Default: "
           << (allele_frequency_usable ? "available" : "not available")
           << "\n";
    }
    for (auto const &tag: _settings.getAlleleFrequencyTags()) {
      tag_type = hdr.getInfoType(tag);
      bool current_tag_allele_frequency_usable = (tag_type == 2);
      bool current_tag_allele_frequency_available = (tag_type != -1);
      allele_frequency_usable |= current_tag_allele_frequency_usable;
      // if (current_tag_allele_frequency_usable) {
      //   usable_tags.insert(tag);
      // }
      if (_settings.warn()) {
        if (current_tag_allele_frequency_available && !current_tag_allele_frequency_usable) {
          cerr << "WARNING: The '" << tag << "' variant allele frequency type ('"
               << ((tag_type == 1) ? "int" : "string")
               << "') is not handled.\n"
               << "         Only tags of type Float are handled.\n";
        }
        cerr << "    - " << tag << ": "
             << (current_tag_allele_frequency_usable ? "available" : "not available")
             << endl;
      }
    }
    if (!allele_frequency_available && _settings.warn()) {
      cerr << "WARNING: No usable variant allele frequency tag found in the file " << f << "\n"
           << "         This is not detrimental, but indentification metric won't take into account variants from this file"
           << endl;
    }

    if (_settings.warn()) {
      cerr << "  - Number of sequences: " << hdr.getSeqnames().size() << endl;
      cerr << "  - Number of samples: " << hdr.nSamples() << " [";
      bool first = true;
      for (auto const &s: hdr.getSamples()) {
        cerr << (first ? "" : ", ") << "'" << s << "'";
      }
      cerr << "]" << endl;
      for (auto const &s: hdr.getSeqnames()) {
        if (!VariantKmerEnumerator::getIndex().getSequenceID(s)) {
          cerr << "WARNING: Sequence '" << s << "' is not indexed.\n"
               << "         Thus k-mer variant index may be incomplete."
               << endl;
        }
      }
    }

    vcfpp::BcfRecord variant(vcf.header); // construct a variant record

    kim_index.unfreeze();
    assert(!kim_index.frozen());
    size_t variants_cpt = 0;
    static const size_t variants_cpt_mask1 = (1 << 4) - 1; // arbitrarly, each 16 added variant we'll update the progress counter
    static const size_t variants_cpt_mask2 = (1 << 20) - 1; // arbitrarly, each 1048576 added variant we'll compress the graph
    while (vcf.getNextVariant(variant)) {
#ifdef DEBUG
      cerr << "Variant: '" << variant << "'" << endl;
#endif
      if (_variantOfInterest(variant)) {
#ifdef DEBUG
        cerr << "Processing variant " << variant.ID()
             << " at pos " << variant.POS()
             << " of sequence " << variant.CHROM()
             << " having length " << (variant.End() - variant.Start())
             << ": " << variant.REF() << " => " << variant.ALT() << endl;
#endif
        VariantKmerEnumerator vke(variant);

        // Compute allele frequencies for all tags
        vaf.compute(vke);
        _dump2temporaryVafFiles(vaf);

        while (vke.nextVariantKmer()) {
          KmerVariantGraph::Edge e = vke.getCurrentKmerVariantEdge();
          bool skip = false;
          vector<pair<string, regex>>::const_iterator it = _kmer_filters.begin();
          while (!skip && (it != _kmer_filters.end())) {
            skip = regex_search(e.kmer, it->second);
            if (!skip) {
              ++it;
            }
          }
          if (skip) {
            if (_settings.warn()) {
              cerr << "    - "
                   << "Skipping k-mer '" << e.kmer << "'"
                   << " since it matches the /" << it->first << "/"
                   << " exclude regular expression." << endl;
            }
          } else {
#ifdef DEBUG
            cerr << "Adding edge: " << e.kmer << " -" << e.rank << "-> " << e.variant << endl;
#endif
            kim_index += e;
          }
        }
      }
      if (_settings.warn()) {
        ++variants_cpt;
        if (!(variants_cpt & variants_cpt_mask1)) {
          cerr << "\033[s\033[B"
               << "  - Number of indexed variants: " << variants_cpt
               << "\033[u" << flush;
        }
        if (!(variants_cpt & variants_cpt_mask2)) {
          cerr << "\033[s\033B"
               << "  - Compressing the index"
               << "\033[u" << flush;
          kim_index.freeze();
          kim_index.unfreeze();
        }
      }
    }
    if (_settings.warn()) {
      cerr << "  - Number of indexed variants: " << variants_cpt << endl;
      cerr << "  - Compressing the index" << endl;
    }
    kim_index.freeze();
    assert(kim_index.frozen());
  }
  assert(kim_index.frozen());
  _closeVafFiles();

  kim_index.unfreeze();
  assert(!kim_index.frozen());
  for (auto const &f: _dna_files) {
    size_t nb_kmers = 0;
    if (_settings.warn()) {
      const DNAFileIndex &index = VariantKmerEnumerator::getIndex();
      nb_kmers = index.getNumberOfNucleotides(f);
      cerr << "* Scanning file '" << f << "'"
           << " (having " << fixed << nb_kmers << " nucleotides)"
           <<" for existing " << _settings.k() << "-mers associated to variant(s)." << endl;
      nb_kmers -= index.getNumberOfSequences(f) * (_settings.k() - 1);
    }
    DNAFileReader reader(_settings.k(), f, false);
    size_t kmer_cpt = 0;
    static const size_t kmer_cpt_mask = (1 << 20) - 1;
    for (string kmer = reader.getNextKmer(true /* skip degenerate */); !kmer.empty(); kmer = reader.getNextKmer(true /* skip degenerate */)) {
      kim_index.setInReferenceKmer(kmer);
      if (_settings.warn() && !(++kmer_cpt & kmer_cpt_mask)) {
        cerr << "\033[s\033[B"
             << "  - Number of processed k-mers: " << kmer_cpt << "/" << nb_kmers
             << "\033[u" << flush;
      }
    }
    if (_settings.warn()) {
      cerr << "\033[s\033[B"
           << "  - Number of processed k-mers: " << nb_kmers
           << endl;
    }
  }
  monitor.stop();
  string metadata;
  metadata += "- Variant files:\n";
  for (auto const &f: _variants_files) {
    metadata += "  - ";
    metadata += f;
    metadata += '\n';
  }
  if (!_variant_filters.empty()) {
    metadata += "- Variant filter(s):\n";
    for (auto const &f: _variant_filters) {
      metadata += "  - ";
      metadata += f;
      metadata += '\n';
    }
  }
  if (!_kmer_filters.empty()) {
    metadata += "- k-mer exclude filter(s):\n";
    for (auto const &f: _kmer_filters) {
      metadata += "  - /";
      metadata += f.first;
      metadata += "/\n";
    }
  }
  metadata += "- Reference files:\n";
  for (auto const &f: _dna_files) {
    metadata += "  - ";
    metadata += f;
    metadata += '\n';
  }
  metadata += "- Resources used for index creation:\n";
  metadata += "  - Wallclock time (in seconds): ";
  metadata += to_string(monitor.getWallClockTime() / 1s);
  metadata += "\n";
  metadata += "  - User CPU time (in seconds): ";
  metadata += to_string(monitor.getUserTime() / 1s);
  metadata += "\n";
  metadata += "  - System CPU time (in seconds): ";
  metadata += to_string(monitor.getSystemTime() / 1s);
  metadata += "\n";
  metadata += "  - Memory: ";
  metadata += Monitor::memoryWithUnit2string(monitor.getMemory());
  metadata += "\n";
  kim_index.extraMetadata(metadata);

  kim_index.freeze();
  assert(kim_index.frozen());

  if (_settings.warn()) {
    cerr << "* Dumping the index to directory"
         << " '" << _settings.getIndexDirectory() << "'"
         << endl;
  }
  kim_index.dump(_settings.getIndexDirectory(), _settings.allowOverwrite());
  _saveTemporaryVafFiles();
  _deleteTemporaryDirectory();
}

int KimProgram::run() {
  try {
    switch (_running_mode) {
    case QUERY_MODE:
      _runQuery();
      break;
    case INDEX_CREATION_MODE:
      _createIndex();
      break;
    default:
      throw Exception("A Bug occurs. Please contact the authors of this program");
    }
  } catch (const exception &e) {
    cerr << "The program has encountered the following error:" << endl
         << endl
         << " => " << e.what() << endl
         << endl;
    return 1;
  };

  if (_settings.warn()) {
    cerr << "Ending time: " << timestamp() << endl
         << "That's All, Folks!!!" << endl;
  }

  return 0;
}


int main(int argc, char **argv) {

  KimProgram &app = KimProgram::App(argc, argv);

  const int exit_code = app.run();

  return exit_code;
}
