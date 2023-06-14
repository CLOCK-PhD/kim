#ifndef __VARIANT_KMER_ASSOCIATION_H__
#define __VARIANT_KMER_ASSOCIATION_H__


#include <string>

class VariantKmerAssociation {

  private: 
    std::string kmer_suffix;
    std::string rs_id;
    std::string chromosome;
    size_t kmer_rank;
    bool in_genome;

  public:
  
    /**
     * Creates an association between the given k-mer suffix and the given variant ID.
     *
     * \param kmer_suffix Suffix of the k-mer (prefixes are not stored since associations are mapped to the given prefix).
     *
     * \param rs_id Variant index.
     *
     * \param kmer_rank Rank of the k-mer for the given variant (from 5' to 3').
     *
     * \param in_genome True if the k-mer is natively present in the genome.
     */
    VariantKmerAssociation(std::string kmer_suffix, std::string rs_id, size_t kmer_rank, bool in_genome);

    /**
     * Get the k-mer suffix.
     *
     * \return Returns the k-mer suffix of this association.
     */
    const std::string &getKmerSuffix() const;

    // Todo : all the getters.    

};

#endif
