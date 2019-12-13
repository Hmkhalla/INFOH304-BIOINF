#define NOMINMAX 
#ifndef DEF_ALGORITHM
#define DEF_ALGORITHM
#include "Database.h"
#include "Protein.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;

class Algorithm 
{
	Database *db;
	Protein *pr;
    int8_t m_blosumMatrix[28][28];
    uint16_t m_extGap;
    uint16_t m_openGap;
    //uint16_t max(uint16_t n1, uint16_t n2);


public:
    Algorithm(string dbPath, string fasta, string pathBlosum);
    ~Algorithm();
    uint16_t max(uint16_t n1, uint16_t n2);
    int8_t getBlosumMatrix();
    //void exactMatch();
    int scoring(uint8_t *seq1, uint32_t nb1, uint8_t *seq2, uint32_t nb2); 
    void notExactMatch();
};
#endif
