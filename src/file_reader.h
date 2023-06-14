#ifndef __FILE_READER_H__
#define __FILE_READER_H__

#include <cstdlib>
#include <string>
#include <fstream>

class FileReader: public std::ifstream {

  private:
  
    const char *filename;
    std::ifstream ifs;
    size_t line, col;
    size_t k;
    std::string read_id;
    std::string kmer;
    size_t kmer_pos;

  public:

    /**
     * Creates a (fastq) file reader.
     *
     * The file is open as an input file stream.
     *
     * \param filename The name of the fastq file to read.
     *
     * \param k The length of the k-mers to extract.
     */
    FileReader(const char *filename, size_t k);
    
    /**
     * Closes the file stream before destruction.
     */
    ~FileReader();

    /**
     * Get the next k-mer from this file reader object.
     *
     * This method may update the current read id and updates the current k-mer position.
     *
     * \return Returns the next available k-mer or an empty string.
     */    
    const std::string &getNextKmer();

    /**
     * Get the current k-mer extracted by this file reader object.
     *
     * \return Returns the current available k-mer or an empty string.
     */    
    const std::string &getCurrentKmer() const;


    /**
     * Get the read ID the current k-mer comes from.
     *
     * \return Returns the current available read ID or an empty string.
     */    
    const std::string &getCurrentRead() const;

    /**
     * Get the position of the current k-mer in the sequence (starting from 0).
     *
     * \return Returns the position of the current available k-mer or -1.
     */    
    size_t getCurrentKmerRelativePosition() const;

};

#endif