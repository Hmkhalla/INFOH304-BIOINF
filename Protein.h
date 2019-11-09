#ifndef DEF_PROTEIN
#define DEF_PROTEIN

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

class Protein {
	private:
	string m_proteinDescription;
	//type? m_proteinSequence;
	public:
	Protein();
	string getProteinDescription();
	//type? getProteinSequence();
};
#endif
