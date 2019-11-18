#ifndef DEF_DATABASE
#define DEF_DATABASE

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

class Database {
	string path;
	uint32_t version;
	uint32_t db_type;
	char * title;
	char * time;
	uint32_t nb_seq;
	uint32_t max_seq;
	uint64_t res_count;
	uint32_t * seq_table;
	uint32_t * head_table;
	public:
	Database(char* pinFile);
	~Database();
	void printDbDescription() const;
	const uint32_t* getSeqTable() const;
	const uint32_t* getHeadTable() const;
	uint32_t getIndexSeq(int index) const;
	uint32_t getIndexHead(int index) const;
	uint8_t* find_seq(int index) const;
	string find_header(int index) const;
	
};
#endif
