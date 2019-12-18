#ifndef DEF_ALGORITHM
#define DEF_ALGORITHM
#include "Database.h"
#include "Sequence.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <thread>

#define BLOSUM_SIZE 24
using namespace std;


class Algorithm 
{
	Database *db;
	int Q;
	int R;
	uint nb_seq;
	uint8_t nb_thread;
    int8_t blosumMatrix[28][28] ={0};
    Sequence query;
    thread* threads;
    Sequence *seqArray;
    map<char, uint8_t> conversion;
    
    int Max(int n1, int n2, int n3) const;
    void initBlosumMatrix(string &blosumPath);
    int scoring(Sequence& tmp_seq) const;


public:
    Algorithm(string dbPath, string fasta, string pathBlosum, int core, int extension_gap = 1, int open_gap = 11);
    ~Algorithm();
    void startMultithread();
    void exactMatch();
    void swAlgo(int start);
    void showResult();
};
#endif
