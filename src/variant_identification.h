#ifndef __VARIANT_IDENITIFICATION_H__
#define __VARIANT_IDENITIFICATION_H__


#include <string>
#include <unordered_map>
#include <list>

class VariantIdentification {

  public:

    /**
     * The type of the read identifier.
     */
    typedef std::string ReadID_type;

  private:

    std::unordered_multimap<ReadID_type, size_t> informations;

  public:

    /**
     * Builds a variant identification object with no informations.
     */
    VariantIdentification();

    /**
     * Add some information associated to some variant.
     *
     * This method consider that a new k-mer (at the given position) is associated to some potential variant for the specified read ID.
     *
     * \param read_id The read header ID.
     *
     * \param pos The position of the k-mer potentially associated to some variation.
     *
     * \return Return the updated variant information.
     */
    VariantIdentification &add(ReadID_type read_id, size_t pos);

    /**
     * Get the list of read IDs.
     *
     * \return Return the list of the read IDs.
     */
    std::list<ReadID_type> getReads() const;

    /**
     * Get the score associated to the given read for the variant detection.
     *
     * TODO: explain how the score is computed and whqat it means.
     *
     * \param read_id The read ID.
     *
     * \return Returns the computed score value. 
     */
    double getReadScore(ReadID_type read_id) const;

};

#endif
