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

#include <iostream>
#include <map>
#include <regex>

#ifdef _OPENMP
#  include <omp.h>
#endif

using namespace std;
using namespace kim;


struct Arg: public option::Arg {

  static void printError(const char* msg1, const option::Option& opt, const char* msg2);

  static option::ArgStatus Unknown(const option::Option& option, bool msg);

  static option::ArgStatus Required(const option::Option& option, bool msg);

  static option::ArgStatus NonEmpty(const option::Option& option, bool msg);

  static option::ArgStatus Numeric(const option::Option& option, bool msg);

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


class KimProgram {

private:

  static const option::Descriptor _usage[];

  struct _OptionHandler {
    option::Stats  stats;
    option::Option *options;
    option::Option *buffer;
    option::Parser parse;

    _OptionHandler(int argc, char **argv);
    ~_OptionHandler();
  };

  enum  _OptionIndex {
    UNKNOWN_OPT,
    HELP_OPT,
    COPYRIGHT_OPT,
    VERSION_OPT,
    QUIET_OPT,
    FORCE_OPT,
    INDEX_DIRECTORY_OPT,
    CREATE_INDEX_OPT,
    VARIANT_FILTER_OPT,
    KMER_LENGTH_OPT,
    KMER_PREFIX_LENGTH_OPT,
    REFERENCE_OPT,
    VARIANTS_OPT,
    CONSISTENCY_CHECKING_OPT,
    OUTPUT_OPT
  };

  enum _RunningMode {
    UNDEFINED_MODE,
    QUERY_MODE,
    INDEX_CREATION_MODE
  };

  string _progname;
  Settings _settings;
  _RunningMode _running_mode;
  vector<string> _variant_filters;
  vector<string> _dna_files;
  vector<string> _variants_files;
  string _output_file;

  KimProgram(int argv, char **argc);
  ~KimProgram();

  void _processOptions(_OptionHandler &_opts);
  bool _variantOfInterest(vcfpp::BcfRecord &v) const;
  void _runQuery();
  void _createIndex();

  static bool _isURL(const string &f);
  static bool _checkIfFileExists(const string &f);
  static void _assertFileExists(const string &f);
  static void _dumpFile(const string &filename, size_t n = size_t(-1), ostream &os = cout);
  static void _dumpResults(const map<string, VariantIdentification> &results,
                            const KmerVariantGraph &index,
                            ostream &os = cout);
  static void _showCopyright(bool full);

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
#define KIM_DEFAULT_KMER_LENGTH         27
#define KIM_DEFAULT_KMER_PREFIX_LENGTH  6

#define _str(x) #x
#define stringify(x) _str(x)

const option::Descriptor KimProgram::_usage[] =
  {
   { UNKNOWN_OPT,            0, "" , "",                   Arg::None,
     PACKAGE " version " VERSION " -- " PACKAGE_SHORT_DESCRIPTION "\n\n"
     "Usages:\n"
     "  kim [options] --create-index\n"
     "  kim [options] <file> [<file> ...]\n"
     "  kim [information option]\n"
     "\n" },
   ////////////////////////////////////////////////////////////////////////
   { UNKNOWN_OPT,              0, "" , "",                   Arg::None,
     "Available information options:" },
   { HELP_OPT,                 0, "h", "help",               Arg::None,
     "  -h | --help \tPrint usage and exit." },
   { VERSION_OPT,              0, "v", "version",            Arg::None,
     "  -v | --version \tPrint version and exit." },
   { COPYRIGHT_OPT,            0, "" , "copyright",          Arg::None,
     "       --copyright \tPrint copyright summary and exit." },
   { COPYRIGHT_OPT,            1, "" , "full-copyright",     Arg::None,
     "       --full-copyright \tPrint full copyright and exit." },
   ////////////////////////////////////////////////////////////////////////
   { UNKNOWN_OPT,              0, "" , "",                   Arg::None,
     "\n"
     "Options available for both index creation and querying:" },
   { QUIET_OPT,                0, "q", "quiet",              Arg::None,
     "  -q | --quiet \tDon't produce warning or extra informations on"
     " standard error channel." },
   { FORCE_OPT,                0, "f", "force",              Arg::None,
     "  -f | --force \tForce overwriting existing index directory (on index"
     " creation) or output file (on querying). This is dangerous, you are"
     " advertised." },
   { INDEX_DIRECTORY_OPT,      0, "d", "index-dir",          Arg::Required,
     "  -d | --index-dir <dir> \tDirectory containing the index files"
     " (default: " KIM_DEFAULT_INDEX_DIRECTORY ")." },
   { CONSISTENCY_CHECKING_OPT, 0, "", "check-consistency",   Arg::None,
     "       --check-consistency \tEnable consistency checking while parsing"
     " DNA sequence files. This significantly increases the program running"
     " time but ensures that your file are correctly written (useful when"
     " something wrong occurs, disabled by default)" },
   { UNKNOWN_OPT,              0, "" , "",                   Arg::None,
     "\n"
     "Options available only for index creation:" },
   { CREATE_INDEX_OPT,         0, "c", "create-index",       Arg::None,
     "  -c | --create-index \tFlag to run in index creation mode. If this"
     " flag is given, the --reference and --variants options must be"
     " provided." },
   { KMER_LENGTH_OPT,          0, "k", "kmer-length",        Arg::Numeric,
     "  -k | --kmer-length <length> \tLength of the k-mers"
     " (default: " stringify(KIM_DEFAULT_KMER_LENGTH) ")." },
   { KMER_PREFIX_LENGTH_OPT,   0, "p", "kmer-prefix-length", Arg::Numeric,
     "  -p | --kmer-prefix-length <length> \tLength of the k-mers prefix"
     " (default: " stringify(KIM_DEFAULT_KMER_PREFIX_LENGTH) ")." },
   { REFERENCE_OPT,            0, "r", "reference",          Arg::Required,
     "  -r | --reference <file> \tReference sequence from which k-mers are"
     " built. It is possible to use several references and at least one"
     " reference file must be provided (see option --variants). The given"
     " reference file must be fasta or fastq formatted." },
   { VARIANTS_OPT,             0, "v", "variants",           Arg::Required,
     "  -v | --variants <file> \tVariants to index. The sequence name of"
     " each variant must be defined in exactly one of the provided reference"
     " files (see option --reference). The variants file must be VCF of BCF"
     " formatted (compressed with gzip or not). You also can provide URL "
     "instead of some local file name." },
   { VARIANT_FILTER_OPT,       0, "F", "variant-filter",     Arg::Required,
     "  -F | --variant-filter <expr> \tFilter variants according to the"
     " given expression. When this option is provided multiple time, each"
     " expression must be satisfyed for the variant to be added into the"
     " index. See filter syntax explanation below for additional details." },
   ////////////////////////////////////////////////////////////////////////
   { UNKNOWN_OPT,              0, "" , "",                   Arg::None,
     "\n"
     "Options available only for Index querying:" },
   { OUTPUT_OPT,               0, "o", "output",             Arg::Required,
     "-o | --output <file> \tPath to the filename where results are stored"
     " (instead of being prompted to the standard output)." },
   ////////////////////////////////////////////////////////////////////////
   { UNKNOWN_OPT,              0, "" ,  "",                  Arg::None,
     "\n"
     PACKAGE_FULL_DESCRIPTION "\n"
     "The kim program has two major features:\n"
     "- The first one is to create a k-mer index related to some variant of"
     " interest.\n"
     "- The second one is to analyse some given input files (containing raw"
     " biological data from some experiment) and to query an existing index"
     " to exhibit variants that are present in the input data.\n"
     "\n"
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
     "INDEX QUERY\n"
     "===========\n"
     "\n"
     "To query an existing index, the kim program needs an index and some"
     " files to analyse. An example of command line to query some index is:\n"
     "\n  kim --index-dir /path/to/my/index/ file1.fastq file2.fastq\n"
     "\n"
     "As for the index creation, the biological sequence files must be either"
     " fasta or fastq formatted."
     "\n" },
     {0, 0, 0, 0, 0, 0}
  };

KimProgram::_OptionHandler::_OptionHandler(int argc, char **argv):
  stats(true, _usage, argc, argv),
  options(new option::Option[stats.options_max]),
  buffer(new option::Option[stats.buffer_max]),
  parse(_usage, argc, argv, options, buffer) {
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
  _variant_filters(),
  _dna_files(), _variants_files(), _output_file()
{
  _settings.setIndexDirectory(KIM_DEFAULT_INDEX_DIRECTORY);

  argc -= (argc > 0);
  argv += (argc > 0); // skip program name argv[0] if present

  if (argc == 0) {
    option::printUsage(cout, _usage);
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

void KimProgram::_processOptions(_OptionHandler &_opts) {

  if (_opts.parse.error()) {
    throw Exception("Parse error");
  }

  if (_opts.options[HELP_OPT]) {
    option::printUsage(cout, _usage);
    exit(0);
  }

  if (_opts.options[VERSION_OPT]) {
    cout << "kim version " << VERSION << endl;
    exit(0);
  }

  if (_opts.options[COPYRIGHT_OPT]) {
    _showCopyright(_opts.options[COPYRIGHT_OPT].type());
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
    _settings.setIndexDirectory(_opts.options[INDEX_DIRECTORY_OPT].arg,
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

    // Check if all given variant files are readable and fill the
    // variant_files vector with their filenames.
    _variant_filters.reserve(_opts.options[VARIANT_FILTER_OPT].count());
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
      _assertFileExists(opt->arg);
      _dna_files.push_back(opt->arg);
    }

    // Check if all given variant files are readable and fill the
    // variant_files vector with their filenames.
    _variants_files.reserve(_opts.options[VARIANTS_OPT].count());
    for (option::Option* opt = _opts.options[VARIANTS_OPT]; opt; opt = opt->next()) {
      if (!_isURL(opt->arg)) _assertFileExists(opt->arg);
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
    if (_opts.options[OUTPUT_OPT]) {
      Exception e;
      e << "Option '" << _opts.options[OUTPUT_OPT].name << "'"
        " is not available for index creation.";
      throw e;
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
      _assertFileExists(_opts.parse.nonOption(i));
      _dna_files.push_back(_opts.parse.nonOption(i));
    }

    // Ensure that given output filename (if provided) is writable.
    if (_opts.options[OUTPUT_OPT]) {
      _output_file = _opts.options[OUTPUT_OPT].arg;
      if (!_settings.allowOverwrite()) {
        if (_checkIfFileExists(_output_file)) {
          Exception e;
          e << "File '" << _output_file << "'"
            << " already exists.\n"
            << "Please use another file name or use option '--force'.";
          throw e;
        }
      }

      ofstream ofs(_output_file);
      if (ofs) {
        ofs.close();
      } else {
        Exception e;
        e << "Unable to open '" << _output_file << "' file to print results.";
        throw e;
      }
    }

    // Prevent using options defined for other mode(s) than index
    // querying.
    if (_opts.options[KMER_LENGTH_OPT]) {
      Exception e;
      e << "Option '" << _opts.options[KMER_LENGTH_OPT].name << "'"
        " is not available for index querying.";
      throw e;
    }
    if (_opts.options[KMER_PREFIX_LENGTH_OPT]) {
      Exception e;
      e << "Option '" << _opts.options[KMER_PREFIX_LENGTH_OPT].name << "'"
        " is not available for index querying.";
      throw e;
    }
    if (_opts.options[REFERENCE_OPT]) {
      Exception e;
      e << "Option '" << _opts.options[REFERENCE_OPT].name << "'"
        " is not available for index querying.";
      throw e;
    }
    if (_opts.options[VARIANTS_OPT]) {
      Exception e;
      e << "Option '" << _opts.options[VARIANTS_OPT].name << "'"
        " is not available for index querying.";
      throw e;
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

bool KimProgram::_checkIfFileExists(const string &f) {
  bool res = false;
  ifstream ifs(f);
  if (ifs) {
    ifs.close();
    res = true;
  }
  return res;
}

void KimProgram::_assertFileExists(const string &f) {
  if (!_checkIfFileExists(f)) {
    Exception e;
    e << "Unable to find or open the file '" << f << "'.";
    throw e;
  }
}

// Dump the file content (up to the number of given characters) to the
// given stream.
void KimProgram::_dumpFile(const string &filename, size_t n, ostream &os) {
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

void KimProgram::_dumpResults(const map<string, VariantIdentification> &results,
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

void KimProgram::_runQuery() {
  Monitor monitor;
  // Loading the index.
  KmerVariantGraph kim_index(_settings);

  kim_index.freeze();
  map<string, VariantIdentification> variants_map;
  // Process each input file
  for (auto const &f: _dna_files) {
    if (_settings.warn()) {
      cerr << "Processing file '" << f << "'" << endl;
    }
    DNAFileReader reader(_settings.k(), f, _settings.warn());

    // Process each k-mer from the current input file
    for (string kmer = reader.getNextKmer(true /* skip degenerate */); !kmer.empty(); kmer = reader.getNextKmer(true /* skip degenerate */)) {
      // cerr << "K-mer : '" << kmer << "'" << endl;
      list<KmerVariantEdgesSubindex::KmerVariantAssociation> assoc_variant = kim_index.search(kmer);
      // cerr << "assoc_variant size is " << assoc_variant.size() << endl;
      for (list<KmerVariantEdgesSubindex::KmerVariantAssociation>::const_iterator it = assoc_variant.cbegin(); it != assoc_variant.cend(); ++it) {
        // cerr << "Current variant association concerns variant '" << it->variant_node.variant << "' which concerns " << it->variant_node.in_degree << " k-mers" << endl;
        // cerr << "This variant is potentially seen in read " << reader.getCurrentSequenceDescription()
        //      << " by the k-mer at position " << reader.getCurrentKmerRelativePosition()
        //      << " in the read" << endl;
        VariantIdentification &v_ident = variants_map[it->variant_node.variant];
        v_ident.add(reader.getCurrentSequenceDescription(), reader.getCurrentKmerRelativePosition());
      }
    }
  }

  // Write the results
  if (!_output_file.empty()) {
    ofstream ofs(_output_file);
    if (ofs) {
      _dumpResults(variants_map, kim_index, ofs);
      ofs.close();
    } else {
      Exception e;
      e << "Unable to open '" << _output_file << "' file to print results.";
      throw e;
    }
  } else {
    _dumpResults(variants_map, kim_index);
  }
  monitor.stop();
  cerr << "# Resources used for this query:" << endl
       << "# - Wallclock time (in seconds): "
       << (monitor.getWallClockTime() / 1s) << endl
       << "# - User CPU time (in seconds): " << (monitor.getUserTime() / 1s) << endl
       << "# - System CPU time (in seconds): " << (monitor.getSystemTime() / 1s) << endl
       << "# - Memory: " << Monitor::memoryWithUnit2string(monitor.getMemory()) << endl;
}

void KimProgram::_createIndex() {

  Monitor monitor;

  Settings fake_settings(_settings);
  try {
    fake_settings.unfreeze();
    fake_settings.setIndexDirectory(_settings.getIndexDirectory(), false, true);
  } catch (BadSettingsException &e) {
    // The index directory already exists.
    if (_settings.allowOverwrite()) {
      fake_settings.setIndexDirectory(_settings.getIndexDirectory()
                                      + "this kind of path should not exist",
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

  VariantKmerEnumerator::init(_settings, _dna_files);

  KmerVariantGraph kim_index(_settings);
  kim_index.freeze();
  assert(kim_index.frozen());
  assert(kim_index.getNbKmerVariantEdges() == 0);
  assert(kim_index.getNbKmers() == 0);
  assert(kim_index.getNbVariants() == 0);

  // Processing Variants files
  for (auto const &f: _variants_files) {
    if (_settings.warn()) {
      cerr << "* Processing variant file '" << f << "'" << endl;
    }
    vcfpp::BcfReader vcf(f);
#ifdef _OPENMP
    vcf.setThreads(omp_get_num_threads());
#endif
    const vcfpp::BcfHeader &hdr = vcf.getHeader();
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
          cerr << "  WARNING: Sequence '" << s << "' is not indexed."
               << " Thus k-mer variant index may be incomplete."
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
        while (vke.nextVariantKmer()) {
          KmerVariantGraph::Edge e = vke.getCurrentKmerVariantEdge();
#ifdef DEBUG
          cerr << "Adding edge: " << e.kmer << " -" << e.rank << "-> " << e.variant << endl;
#endif
          kim_index += e;
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
