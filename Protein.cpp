#include <iostream>
#include <fstream>
#include <string>
#include "Protein.h"
#include <vector>

using namespace std;

    string m_proteinDescription;
    string m_proteinSequence;
    std::vector <int> m_vectorSequence;
    
    Protein:: Protein() {
	ifstream queryFlux("Add path to Fasta file");
	if(queryFlux.is_open()){
	    string line;
	    while (getline(queryFlux,line).good()) {
		    if(line[0]=='>'){
			    m_proteinDescription=line.substr(1);}
		    else{
			    m_proteinSequence+=line;}	
	    }
	}
	else {
		cout <<"Error: Impossible to open the file." <<endl;}
	}


    string Protein::getProteinDescription() {   
	    return m_proteinDescription;}
	    
    string Protein::getProteinSequence() {
		return m_proteinSequence;}

    void Protein::convertToValue(){
		int length=m_proteinSequence.length();
		for( int i=0; i<length; i++) {
			char val=m_proteinSequence[i];
			switch (val) {
				case '-':
				    m_vectorSequence.push_back(0);
				    break;
				case 'A':
				    m_vectorSequence.push_back(1);
				    break;
				case 'B':
				    m_vectorSequence.push_back(2);
				    break;
				case 'C':
				    m_vectorSequence.push_back(3); 
				    break;
				case 'D':
				    m_vectorSequence.push_back(4);
				    break;
				case 'E':
				    m_vectorSequence.push_back(5);
				    break;
				case 'F':
				    m_vectorSequence.push_back(6);
				    break;
				case 'G':
				    m_vectorSequence.push_back(7);
				    break;
				case 'H':
				    m_vectorSequence.push_back(8);
				    break;
				case 'I':
				    m_vectorSequence.push_back(9);
				    break;
				case 'J':
				    m_vectorSequence.push_back(27);
				    break;
				case 'K':
				    m_vectorSequence.push_back(10);
				    break;
				case 'L':
				    m_vectorSequence.push_back(11);
				    break;
				case 'M':
				    m_vectorSequence.push_back(12);
				    break;
				case 'N':
				    m_vectorSequence.push_back(13);
				    break;
				case 'O':
				    m_vectorSequence.push_back(26);
				    break;
				case 'P':
				    m_vectorSequence.push_back(14);
				    break;
				case 'Q':
				    m_vectorSequence.push_back(15);
				    break;
				case 'R':
				    m_vectorSequence.push_back(16);
				    break;
				case 'S':
				    m_vectorSequence.push_back(17);
				    break;
				case 'T':
				    m_vectorSequence.push_back(18);
				    break;
				case 'U':
				    m_vectorSequence.push_back(24);
				    break;
				case 'V':
				    m_vectorSequence.push_back(19);
				    break;
				case 'W':
				    m_vectorSequence.push_back(20);
				    break;
				case 'X':
				    m_vectorSequence.push_back(21);
				    break;
				case 'Y':
				    m_vectorSequence.push_back(22);
				    break;
				case 'Z':
				    m_vectorSequence.push_back(23);
				    break;
				case '*':
				    m_vectorSequence.push_back(25);
				    break;
				default:
				    cout<<"Error"<<endl;
				    break;
		        } 
		    }
		}		 
		
		
    void Protein::printVectorSequence() {
		for (std::vector<int>::const_iterator i= m_vectorSequence.begin(); i!= m_vectorSequence.end(); ++i) {
            cout << *i << "-";}
		}