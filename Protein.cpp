#include <iostream>
#include <fstream>
#include <string>
#include "Protein.h"


using namespace std;

    string m_proteinDescription;
    //type? m_proteinSequence;
    
    Protein:: Protein() {
	ifstream queryFlux("Enter path to FASTA file");
	if(queryFlux.is_open()) {
		getline(queryFlux, m_proteinDescription);
		queryFlux.close();}
	else {
		cout <<"Error: Impossible to open the file." <<endl;}
	}


    string Protein::getProteinDescription() {   
	    return m_proteinDescription;}
    //type? Protein::getProteinSequence() 
