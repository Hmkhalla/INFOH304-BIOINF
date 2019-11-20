#include<iostream>
#include<fstream>
#include "Database.h"
using namespace std;


Database::Database(char* db_fasta_file){
	path = db_fasta_file;
	ifstream in (path+".pin", ios::binary);
	if( in.is_open() )
	{
		
		uint32_t length;
		
		in.read((char *) (&version), sizeof(version) );version=__bswap_32(version);
		in.read((char *) (&db_type), sizeof(db_type) );db_type=__bswap_32(db_type);

		in.read((char *) (&length), sizeof(length) );length = __bswap_32(length);
		title = new char[length+1];
		for (int i=0; i< length; i++) in.read((char *) (&title[i]), sizeof(char)); title[length]=0;
		
		in.read((char *) (&length), sizeof(length) );length=__bswap_32(length);
		time = new char[length+1];
		for (int i=0; i< length; i++) in.read((char *) (&time[i]), sizeof(char)); time[length]=0;
		
		in.read((char *) (&nb_seq), sizeof(nb_seq) );nb_seq=__bswap_32(nb_seq);
		in.read((char *) (&res_count), sizeof(res_count) );
		in.read((char *) (&max_seq), sizeof(max_seq) );max_seq=__bswap_32(max_seq);
		
		Index_seq_table = new uint32_t[nb_seq+1];
		Index_head_table = new uint32_t[nb_seq+1];
		
		for (int i=0; i<= nb_seq; i++){
			in.read((char *) (&Index_head_table[i]), sizeof(int32_t));
			Index_head_table[i]=__bswap_32(Index_head_table[i]);}
		
		for (int i=0; i<= nb_seq; i++){
			in.read((char *) (&Index_seq_table[i]), sizeof(int32_t));
			Index_seq_table[i]=__bswap_32(Index_seq_table[i]);}
		
		in.close();
		prot_table.reserve(nb_seq);
		ifstream psq_f (path+".psq", ios::binary);
		ifstream phr_f (path+".phr", ios::binary);
		for (int i=0; i<nb_seq; i++){
			string m_proteinDescription = find_header((int)getIndexHead(i), (int) getIndexHead(i+1), phr_f);
			prot_table.push_back(Protein(m_proteinDescription, find_seq(i, psq_f)));
		}
		psq_f.close();
		phr_f.close();
	}
	else
		cerr << "Impossible d'ouvrir le fichier pin" << endl;
}

Database::~Database(){
	delete title;
	delete time;
	delete[] Index_seq_table;
	delete[] Index_head_table;
	prot_table.clear();
}

void Database::printDbDescription() const{
	cout << "Version :" << version << endl;
	cout << "db_type :" << db_type << endl;
	cout << "title :" << title << endl;
	cout << "Timestamp :" << time << endl;
	cout << "Number of sequences :" << nb_seq << endl;
	cout << "Residue count :" << res_count << endl;
	cout << "Maximum sequence :" << max_seq << endl;
	for (int i =0; i< 5; i++){
		cout << "Header offset table " << i << " : "<< Index_head_table[i] << endl;
	}
	cout << "Header offset table " << nb_seq << " : "<< Index_head_table[nb_seq] << endl;
	cout << "Sequence offset table " << nb_seq << " : "<< Index_seq_table[nb_seq] << endl;
	for (int i =0; i< 5; i++){
		cout << "Sequence offset table " << i << " : "<< Index_seq_table[i] << endl;
	}
}

const uint32_t* Database::getSeqTable() const{
	return Index_seq_table;
}

const uint32_t* Database::getHeadTable() const{
	return Index_head_table;
}

uint32_t Database::getIndexSeq(int index) const{
	if (index > nb_seq){cerr << "Seq index "<< index<< " out of range"<< endl; return 0;}
	return Index_seq_table[index];
}

uint32_t Database::getIndexHead(int index) const{
	if (index > nb_seq){cerr << "Head index "<< index<< " out of range"<< endl; return 0;}
	return Index_head_table[index];
}


vector <uint8_t> Database::find_seq(int index, ifstream &in) const{
	int len = ((int)getIndexSeq(index+1)-(int)getIndexSeq(index));
	vector<uint8_t> bytes(len);
	if( in.is_open() )
	{
		in.seekg(getIndexSeq(index)*sizeof(uint8_t));
		in.read((char *) &(bytes[0]), sizeof(uint8_t)*len);
	}
	else
		{cerr << "Impossible d'ouvrir le fichier psq" << endl;}
	return bytes;
}


string Database::find_header(int start, int end, ifstream &in) const{
	if( in.is_open() ){
	string title, type, db, id, res;
	bool readingInt =false;
	bool readingLen =false;
	bool readingStr =false;
	uint8_t byte;
	in.seekg(start*sizeof(byte));
	while( in.tellg() < end && in.read((char *) (&byte), sizeof(byte))){
		if (! readingInt && ! readingLen && !readingStr){
			switch(byte) {
			  case 0x1a:
				// String
				readingStr = true;
				readingLen = true;
				break;
			  case 0x02:
				// int
				readingInt = true;
				readingLen = true;
				break;
			}
		}
		else if (readingLen){
			if (readingStr){
				int len = 0;
				if (byte & (1<<7)){
					byte = (byte & ~(1 << 7)) | ((0 << 7) & (1 << 7));
					in.read((char *) &(len), sizeof(byte)*(int)byte);
				}else len = (int) byte;
				char buffer[len+1];
				in.read((char *) &(buffer[0]), sizeof(byte)*len);
				buffer[len] = 0;
				readingStr = false;
				readingLen = false;
				if (title.empty()) title = buffer;
				else db = buffer;
			}else if (readingInt){
				int len = (int) byte;
				uint buffer;
				in.read((char *) &(buffer), sizeof(byte)*len);
				readingInt = false;
				readingLen = false;
				if (id.empty()) id = to_string((int)buffer);
				else type = to_string((int)buffer);
			}
		}		
	}
	res = type +"|"+db + "|" + id + " " + title;
	return res;
	}
	else cerr << "Impossible d'ouvrir le fichier" << endl;
}
