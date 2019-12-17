#ifndef DEF_DATABASE
#define DEF_DATABASE

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
using namespace std;

class Database
{
    uint8_t * psq;
    uint32_t *Index_seq_table;
    
    uint8_t * phr;
    uint32_t *Index_head_table;
    
    uint32_t nb_seq;
    
    string description;
    
    uint32_t getIndexSeq(int index) const;
    uint32_t getIndexHead(int index) const;
    void readPin(string & path);
    void loadPhr(string & path);
    void loadPsq(string & path);

public:
    Database(string &dbPath);
    ~Database();
    void printDbDescription() const;
    
    uint8_t *getSeq(int index, int &nb) const;
    uint8_t *getHeader(int index, int &nb) const;
    void find_header(string &res, int index) ;
    
    uint32_t getNbSeq() const;
    uint32_t getLenSeq(int index) const;
};
#endif
