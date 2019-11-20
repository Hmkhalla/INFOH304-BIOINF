#ifndef DEF_PROTEIN
#define DEF_PROTEIN

#include <iostream>
#include <fstream>
#include <string>
#include <vector>


using namespace std;

class Protein {
	private:
	string m_proteinDescription;
	std::vector<uint8_t> m_vectorSequence;
	
	public:
	Protein();
	~Protein();
	Protein(string m_proteinDescription, vector<uint8_t> m_vectorSequence);
	string getProteinDescription();
	string getProteinSequence();
	void convertToValue();
	void printVectorSequence();
};
#endif
