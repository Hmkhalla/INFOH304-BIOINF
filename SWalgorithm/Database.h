#ifndef DEF_DATABASE
#define DEF_DATABASE

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
//#include <algorithm>
#include <map>
using namespace std;

class Database
{
    ifstream phr;
    ifstream pin;
    ifstream psq;
    uint32_t version;
    uint32_t db_type;
    uint32_t nb_seq;
    uint32_t max_seq;
    uint64_t res_count;
    char *title;
    char *time;
    uint32_t *Index_seq_table;
    uint32_t *Index_head_table;

public:
    Database(string dbPath);
    ~Database();
    void printDbDescription() const;
    uint32_t getIndexSeq(int index) const;
    uint32_t getIndexHead(int index) const;
    uint8_t *find_seq(int index, int &nb);
    void find_header(string &res, int index);
    uint32_t getNbSeq();
};
#endif
