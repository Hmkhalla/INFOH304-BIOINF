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
		
		seq_table = new uint32_t[nb_seq+1];
		head_table = new uint32_t[nb_seq+1];
		
		for (int i=0; i<= nb_seq; i++){
			in.read((char *) (&head_table[i]), sizeof(int32_t));
			head_table[i]=__bswap_32(head_table[i]);}
		
		for (int i=0; i<= nb_seq; i++){
			in.read((char *) (&seq_table[i]), sizeof(int32_t));
			seq_table[i]=__bswap_32(seq_table[i]);}
		
		in.close();
	}
	else
		cout << "Impossible d'ouvrir le fichier pin" << endl;
}

Database::~Database(){
	delete title;
	delete time;
	delete[] seq_table;
	delete[] head_table;
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
		cout << "Header offset table " << i << " : "<< head_table[i] << endl;
	}
	cout << "Header offset table " << nb_seq << " : "<< head_table[nb_seq] << endl;
	cout << "Sequence offset table " << nb_seq << " : "<< seq_table[nb_seq] << endl;
	for (int i =0; i< 5; i++){
		cout << "Sequence offset table " << i << " : "<< seq_table[i] << endl;
	}
}

const uint32_t* Database::getSeqTable() const{
	return seq_table;
}

const uint32_t* Database::getHeadTable() const{
	return head_table;
}

uint32_t Database::getIndexSeq(int index) const{
	if (index >= max_seq){cout << "index "<< index<< " out of range"<< endl; return 0;}
	return seq_table[index];
}

uint32_t Database::getIndexHead(int index) const{
	if (index >= max_seq){cout << "index "<< index<< " out of range"<< endl; return 0;}
	return head_table[index];
}

uint8_t * Database::find_seq(int index) const{
	if (index >= max_seq){cout << "index "<< index<< " out of range"<< endl; return NULL;}
	ifstream in (path+".psq", ios::binary);
	uint8_t *bytes = NULL;
	if( in.is_open() )
	{
		int len = getIndexSeq(index+1)-getIndexSeq(index);
		bytes= new uint8_t[len];
		in.seekg(getIndexSeq(index)*sizeof(uint8_t));
		for (int i = 0; i<len; i++) in.read((char *) (&bytes[i]), sizeof(uint8_t) );
		in.close();
	}
	else
		cout << "Impossible d'ouvrir le fichier psq" << endl;
	return bytes;
}

//string find_header(int index) const
