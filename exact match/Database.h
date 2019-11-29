#ifndef DEF_DATABASE
#define DEF_DATABASE

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
//#include <algorithm>
using namespace std;

class Database {
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
	Database(string path);
	~Database();
	void printDbDescription() const;
	uint32_t getIndexSeq(int index) const;
	uint32_t getIndexHead(int index) const;
	void find_seq(vector<uint8_t> &bytes, int index);
	void find_header(string &res, int index);
	void exactMatch(char *queryPath);
	//uint8_t scoring(vector<uint8_t> *seq1, uint32_t nb1, vector<uint8_t> *seq2, uint32_t nb2);
	//void notExactMatch(char *queryPath); 
};

void convertToValue(vector<uint8_t> &sequence, ifstream &in);
#endif
