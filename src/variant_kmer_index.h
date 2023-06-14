#ifndef __VARIANT_KMER_INDEX_H__
#define __VARIANT_KMER_INDEX_H__

#include <cstdlib>
#include <list>
#include <map>
#include <string>
#include <variant_identification.h>
#include <variant_kmer_association.h>


class VariantKmerIndex {

  private:
  
    size_t k, k1, k2;
    
    std::map<size_t, VariantKmerAssociation> index;

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
    std::list<VariantIdentification::key_type> search(const std::string &kmer) const;    

};


#endif
