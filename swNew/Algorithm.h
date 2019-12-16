#ifndef DEF_ALGORITHM
#define DEF_ALGORITHM
#include "Database.h"
#include "Sequence.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;

class Algorithm 
{
	Database *db;
	int Q;
	int R;
	uint nb_seq;
    int8_t blosumMatrix[28][28];
    Sequence query;
    map<char, uint8_t> conversion;
    
    int Max(int n1, int n2, int n3) const;
    void initBlosumMatrix(string &blosumPath);
    int scoring(Sequence& tmp_seq) const;


public:
    Algorithm(string dbPath, string fasta, string pathBlosum, int extension_gap = 1, int open_gap = 11);
    ~Algorithm();
    void exactMatch();
    void swAlgo();
};
#endif
