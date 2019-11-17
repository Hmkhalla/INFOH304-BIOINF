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
	string m_proteinSequence;
	std::vector<int> m_vectorSequence;
	
	public:
	Protein();
	string getProteinDescription();
	string getProteinSequence();
	void convertToValue();
	void printVectorSequence();
};
#endif