#ifndef DEF_DATABASE
#define DEF_DATABASE

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
using namespace std;

class Database
{
    uint8_t * phr;
    uint8_t * psq;
    map<char, uint8_t> conversion;
    int8_t blosumMatrix[28][28];
    uint32_t version;
    uint32_t db_type;
    uint32_t nb_seq;
    uint32_t max_seq;
    uint64_t res_count;
    char *title;
    char *time;
    uint32_t *Index_seq_table;
    uint32_t *Index_head_table;
    int16_t Max(int16_t n1, int16_t n2, int16_t n3, int16_t n4);
    void getBlosumMatrix(string blosumPath);

public:
    Database(string path);
    ~Database();
    void printDbDescription() const;
    void readPin(string  path);
    uint32_t getIndexSeq(int index) const;
    uint32_t getIndexHead(int index) const;
    uint8_t *find_seq(int index, int &nb);
    void find_header(string &res, int index);
    void exactMatch(char *queryPath);
    int scoring(uint8_t *seq1, uint32_t nb1, uint8_t *seq2, uint32_t nb2);
    void notExactMatch(char *queryPath);
    uint8_t *convertToValue(ifstream &in, int &nb);
};

#endif
