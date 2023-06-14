#ifndef __VARIANT_KMER_INDEX_H__
#define __VARIANT_KMER_INDEX_H__

#include <cstdlib>
#include <list>
#include <vector>
#include <unordered_map>
#include <string>

#include <variant_identification.h>


class VariantKmerIndex {


  public:

    /**
     * The type of the variant identifier.
     */
    typedef std::string VariantID_type;

    /**
     * The association between some given k-mer and some variant
     */
    struct VariantKmerAssociation {
      /**
       * The variant ID.
       */
      VariantID_type rs_id;
      /**
       * The rank of the k-mer in the variant k-mer markers.
       */
      size_t kmer_rank;
      /**
       * True if the k-mer appears natively in the genome.
       */
      bool in_genome;
    };
    
  private:
  
    size_t k, k1, k2;
    
    std::vector<std::unordered_multimap<VariantID_type, VariantKmerAssociation> > index;

  public:

    /**
     * Load the index of k-mers associated to variants using the files form the given directory.
     *
     * \param path Directory where index files are stored.
     */
    VariantKmerIndex(const char *path);

    /**
     * Get the length of the k-mers.
     *
     * \return Returns the length of the k-mers.
     */
    size_t getKmerLength() const;
    
    /**
     * Get the list of variant identifications involving the given k-mer.
     *
     * \param kmer A string representing the k-mer (k-mer lookup is case insensitive).
     *
     * \return Returns the list of variant identifications involving the given k-mer.
     */
    std::list<VariantKmerAssociation> search(const std::string &kmer) const;

};


#endif
