#include <iostream>
#include <fstream>
#include <string>
#include "Protein.h"
#include <vector>


Protein:: Protein() {
}
Protein:: Protein(string description, vector<uint8_t> sequence) {
	m_proteinDescription = description;
	m_vectorSequence = sequence;
}

Protein::~Protein(){
	m_vectorSequence.clear();
}

string Protein::getProteinDescription() {   
    return m_proteinDescription;}
    
/*
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
	*/	
		
void Protein::printVectorSequence() {
	for (int i = 0; i<m_vectorSequence.size(); i++) {
		cout << (int)m_vectorSequence[i] <<"-";}
		cout << endl;
}
/*	
void Protein::init_seq(vector<uint8_t>& bytes ,int index, ifstream &in) const{
	int len = ((int)getIndexSeq(index+1)-(int)getIndexSeq(index));
	bytes.reserve(len);
	cout << "1" << endl;
	if( in.is_open() )
	{
		in.seekg(getIndexSeq(index)*sizeof(uint8_t));
		cout << "2" << endl;
		in.read((char *) &bytes[0], sizeof(uint8_t)*len);
		cout << "3" << endl;
	}
	else
		{cerr << "Impossible d'ouvrir le fichier psq" << endl;}
	cout << "capa :" << bytes->capacity()<< endl;
	cout << "here" << endl;
	cout << "size :" << bytes->size()<<endl;
	cout << "here:" << (*bytes)[3] <<endl;
}
*/
