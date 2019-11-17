#include <iostream>
#include <fstream>
#include <string>
#include "Protein.h"
#include <vector>

using namespace std;



int main() {
	Protein protein;
	cout<<protein.getProteinDescription()<<endl ;
	cout<<protein.getProteinSequence()<<endl;
	protein.convertToValue();
	protein.printVectorSequence();
	return 0;
};
	