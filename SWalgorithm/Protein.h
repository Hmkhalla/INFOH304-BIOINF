#ifndef DEF_PROTEIN
#define DEF_PROTEIN

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

using namespace std;

class Protein
{
	string m_proteinDescription;
	string m_proteinSequence;
	map<char, uint8_t> m_conversion;
	uint8_t *m_sequenceConverted;

public:
	Protein(string fasta);
	~Protein();
	string getProteinDescription();
	string getProteinSequence();
	uint8_t *convertToValue(string m_proteinSequence);
	uint8_t *getSequenceConverted();
	map<char, uint8_t> getConversion();

};
#endif
