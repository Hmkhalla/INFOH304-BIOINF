#ifndef DEF_PROTEIN
#define DEF_PROTEIN

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

class Protein
{
private:
	string m_proteinDescription;

public:
	std::vector<uint8_t> m_vectorSequence;
	
	public:
	Protein(char* fasta);
	~Protein();
	Protein(string m_proteinDescription, vector<uint8_t> m_vectorSequence);
	Protein(string fasta);
	string getProteinDescription();
	string getProteinSequence();
	void convertToValue();
	void printVectorSequence();
	vector<uint8_t> getVectorSequence();
};
#endif
