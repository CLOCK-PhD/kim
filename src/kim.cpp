#include <iostream>
#include "optionparser.h"
#include "file_reader.h"
#include "variant_kmer_association.h"
#include "variant_kmer_index.h"
#include "variant_identification.h"

using namespace std;


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
  UNKNOWN,
  HELP,
  VERSION,
  INDEX_DIRECTORY,
  INPUT_FILE
};

const option::Descriptor usage[] =
{
  {UNKNOWN,            0, "" , "",          option::Arg::None, "Usage: kim [options] <file> [<file> ...]\n\n"
                                                               "Options:" },
  {HELP,               0, "h", "help",      option::Arg::None, "  -h | --help  \tPrint usage and exit." },
  {VERSION,            0, "v", "version",   option::Arg::None, "  -v | --version  \tPrint version and exit." },
  {INDEX_DIRECTORY,    0, "d", "index-dir", Arg::Required,     "  -d | --index-dir <dir>  \tDirectory containing the index files." },
  {UNKNOWN,            0, "" ,  "",         option::Arg::None, "\nInput files are expected to be fastq formatted"
                                                               "\nExample:\n"
                                                               "  kim --index-dir /path/to/my/index/ file1.fastq file2.fastq\n" },
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

  if (options[HELP] || (argc == 0)) {
    option::printUsage(std::cout, usage);
    return 0;
  }

  if (options[VERSION]) {
    cout << "kim version " << VERSION << endl;
    return 0;
  }

  for (option::Option* opt = options[UNKNOWN]; opt; opt = opt->next()) { 
    std::cout << "Unknown option: " << opt->name << "\n";    
  } 

  if (options[UNKNOWN]) {
    option::printUsage(std::cout, usage);
    return 1;
  }

  const char *index_directory = NULL;
  if (options[INDEX_DIRECTORY]) {
    index_directory = options[INDEX_DIRECTORY].arg;
  } else {
    index_directory = "kim_index";
  }
  
  cout << "Index directory: '" << index_directory << "'" << endl;
  VariantKmerIndex kim_index(index_directory);
  size_t k = kim_index.getKmerLength();
  map<VariantIdentification::key_type, VariantIdentification> variants_map;

  // For each input file
  for (int i = 0; i < parse.nonOptionsCount(); ++i) { 
    cout << "Processing file '" << parse.nonOption(i) << "'" << endl;
    FileReader reader(parse.nonOption(i), k);
    for (std::string kmer = reader.getNextKmer(); !kmer.empty(); kmer = reader.getNextKmer()) {
      list<VariantIdentification::key_type> variant_ids = kim_index.search(kmer);
      for (list<VariantIdentification::key_type>::const_iterator it = variant_ids.cbegin(); it != variant_ids.cend(); ++it) {
        VariantIdentification &v_ident = variants_map[*it];
        v_ident.add(reader.getCurrentRead(), reader.getCurrentKmerRelativePosition());
      }
    }
  }

  cout << "---" << endl
       << "SNPs:" << endl;
  for (map<VariantIdentification::key_type, VariantIdentification>::const_iterator it = variants_map.cbegin();
       it != variants_map.cend();
       ++it) {
    cout << "  - id: " << it->first << endl
         << "    reads:" << endl;
    list<VariantIdentification::read_type> reads = it->second.getReads();
    for (list<VariantIdentification::read_type>::const_iterator read_it = reads.begin();
         read_it != reads.end();
         ++read_it) {
      cout << "      - " << *read_it << ":" << it->second.getReadScore(*read_it) << endl;
    }
  }  

  cout << "That's All, Folks!!!" << endl;
  return 0;
}
