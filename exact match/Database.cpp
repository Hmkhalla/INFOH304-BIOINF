#include "Database.h"
using namespace std;


Database::Database(string path){
	
	pin = ifstream(path+".pin", ios::binary);
	phr = ifstream(path+".phr", ios::binary);
	psq = ifstream(path+".psq", ios::binary);
	if( pin.is_open() )
	{
		
		uint32_t length;
		
		pin.read((char *) (&version), sizeof(version) );version=__bswap_32(version);
		pin.read((char *) (&db_type), sizeof(db_type) );db_type=__bswap_32(db_type);

		pin.read((char *) (&length), sizeof(length) );length = __bswap_32(length);
		title = new char[length+1];
		pin.read((char *) &title[0], sizeof(char)*length); title[length]=0;
		
		pin.read((char *) &length, sizeof(length) );length=__bswap_32(length);
		time = new char[length+1];
		pin.read((char *) &time[0], sizeof(char)*length); time[length]=0;
		
		pin.read((char *) (&nb_seq), sizeof(nb_seq) );nb_seq=__bswap_32(nb_seq);
		pin.read((char *) (&res_count), sizeof(res_count) );
		pin.read((char *) (&max_seq), sizeof(max_seq) );max_seq=__bswap_32(max_seq);
		
		Index_seq_table = new uint32_t[nb_seq+1];
		Index_head_table = new uint32_t[nb_seq+1];
		
		pin.read((char *) &Index_head_table[0], sizeof(int32_t)*(nb_seq+1));
		pin.read((char *) &Index_seq_table[0], sizeof(int32_t)*(nb_seq+1));
		for (int i=0; i<= nb_seq; i++){
			Index_seq_table[i]=__bswap_32(Index_seq_table[i]);
			Index_head_table[i]=__bswap_32(Index_head_table[i]);}
		
		pin.close();
	}
	else
		cerr << "Impossible d'ouvrir le fichier pin" << endl;
}

Database::~Database(){
	pin.close();
	phr.close();
	psq.close();
	delete title;
	delete time;
	delete[] Index_seq_table;
	delete[] Index_head_table;
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


uint32_t Database::getIndexSeq(int index) const{
	if (index > nb_seq){cerr << "Seq index "<< index<< " out of range"<< endl; return 0;}
	return Index_seq_table[index];
}

uint32_t Database::getIndexHead(int index) const{
	if (index > nb_seq){cerr << "Head index "<< index<< " out of range"<< endl; return 0;}
	return Index_head_table[index];
}


void Database::find_seq(vector <uint8_t> &bytes, int index){
	
	int len = ((int)getIndexSeq(index+1)-(int)getIndexSeq(index)) -1; 
	bytes.resize(len);
	if( psq.is_open() )
	{
		psq.seekg((int)getIndexSeq(index)*sizeof(uint8_t), ios::beg);
		psq.read((char *) &(bytes[0]), sizeof(uint8_t)*len);
	}
	else
		{cerr << "Impossible d'ouvrir le fichier psq" << endl;}
}

void Database::find_header(string& res, int index) {
	if( phr.is_open() ){
	string title=""; string type=""; string db=""; string id=""; 
	int start, end; start = getIndexHead(index); end = getIndexHead(index+1);
	bool readingInt =false;
	bool readingLen =false;
	bool readingStr =false;
	uint8_t byte;
	phr.seekg(start*sizeof(byte));
	while( phr.tellg() < end && phr.read((char *) (&byte), sizeof(byte))){
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
			int len = 0;
			if (readingStr){
				if (byte & (1<<7)){
					byte = (byte & ~(1 << 7)) | ((0 << 7) & (1 << 7));
					phr.read((char *) &(len), sizeof(byte)*(int)byte);
				}else len = (int) byte;
				char buffer[len+1];
				phr.read((char *) &(buffer[0]), sizeof(byte)*len);
				buffer[len] = 0;
				readingStr = false;
				if (title.empty()) title = buffer;
				else db = buffer;
			}else if (readingInt){
				len = (int) byte;
				unsigned long intFound =0;
				uint8_t buffer = 0;
				for (int j = 0; j < len; j++){ 
					phr.read((char*)&buffer, sizeof(uint8_t));
					intFound = intFound << 8 | (unsigned long)buffer;
				}
				readingInt = false;
				if (id.empty()) id = to_string((int)intFound);
				else type = to_string((int)intFound);
			}
			readingLen = false;
		}		
	}
	res = type +"|"+db + "|" + id + " " + title;
	}
	else cerr << "Impossible d'ouvrir le fichier" << endl;
}

void Database::exactMatch(char * queryPath){
	vector<uint8_t> query_seq, temp_seq;
	ifstream queryFlux (queryPath, ios::binary);
	convertToValue(query_seq, queryFlux);
	queryFlux.close();
	
	for (int i = 0; i< nb_seq; i++){
		if (getIndexSeq(i+1)-getIndexSeq(i) -1 == query_seq.size()){
			find_seq(temp_seq, i);
			//cout << "admissible :" << i <<endl;
			if (query_seq == temp_seq){
				string header;
				find_header(header, i);
				cout << "answer here: " <<header << endl;
			}
		}
		temp_seq.clear();
	}
}

uint8_t Database::scoring(vector<uint8_t> *seq1, uint32_t nb1, vector<uint8_t> *seq2, uint32_t nb2)
{
	uint8_t gap = -1;
	uint8_t open_gp = -11;
	//vector<int8_t> *v1 = new vector<int8_t>(nb1 + 1, 0);
	vector<uint8_t> *v2 = new vector<uint8_t>(nb2 + 1, 0);
	vector<uint8_t> *curV = new vector<uint8_t>(nb2 + 1, 0);
	uint8_t max_scoring = 0;
	uint8_t Hi_1j_1;
	uint8_t Hij_1;
	uint8_t Hi_1j;
	uint8_t Hij; //plutot les valeurs partantes de ces "cases"
	for (int i = 1; i <= nb1; i++)
	{
		for (int j = 1; i <= nb2 + 1; i++)
		{
			curV->at(j) = max(max(curV->at(j - 1) + open_gp, v2->at(i - 1) + open_gp), max(v2->at(j - 1) + ((seq1->at(i) < seq2->at(i)) ? 1 : 2), 0));
			if (curV->at(j) > max_scoring)
			{
				max_scoring = curV->at(i);
			}
		}
		curV, v2 = v2, curV;
	}
	delete curV;
	delete v2;
	return max_scoring;
}

void Database::notExactMatch(char *queryPath)
{
	vector<uint8_t> query_seq, temp_seq;
	vector<uint8_t> scores;
	uint8_t score;
	vector<uint8_t, uint8_t> classement;
	classement.resize(50);
	ifstream queryFlux(queryPath, ios::binary);
	convertToValue(query_seq, queryFlux);
	queryFlux.close();
	scores.resize(nb_seq);
	for (int i = 0; i < nb_seq; i++)
	{
		find_seq(temp_seq, i);
		score = scoring(&query_seq, nb_seq, &temp_seq, getIndexSeq(i + 1) - getIndexSeq(i) - 1);
		scores.push_back(score);
	}
	uint8_t minMax;
	uint8_t minMaxIndex;
	uint8_t curMax = 0;
	uint8_t curMaxIndex;
	vector<int, int> bestMatches;
	bestMatches.resize(50);
	for (int i = 0; i < 50; i++)
	{
		for (int j = 0; j < nb_seq; j++)
		{
			if (scores.at(j) > curMax && (scores.at(j) < minMax || (scores.at(j) = minMax && j > minMaxIndex)))
			{
				curMax = scores.at(j);
				curMaxIndex = j;
			}
		}
		minMax = curMax;
		minMaxIndex = curMaxIndex;
		bestMatches.at(i) = ((int)curMax, (int)curMaxIndex);
	} // ===> On créer une liste avec les scores de protéines et et leur index dans les offset
}

void convertToValue(vector<uint8_t> &sequence, ifstream &queryFlux)
{
	string buffer;
	if (queryFlux.is_open())
	{
		string line;
		while (getline(queryFlux, line).good())
		{
			if (line[0] != '>')
			{
				buffer += line;
			}
		}
		sequence.reserve(buffer.length() + 1);
		for (int i = 0; i < buffer.length(); i++)
		{
			switch (buffer[i])
			{
			case '-':
				sequence.push_back(0);
				break;
			case 'A':
				sequence.push_back(1);
				break;
			case 'B':
				sequence.push_back(2);
				break;
			case 'C':
				sequence.push_back(3);
				break;
			case 'D':
				sequence.push_back(4);
				break;
			case 'E':
				sequence.push_back(5);
				break;
			case 'F':
				sequence.push_back(6);
				break;
			case 'G':
				sequence.push_back(7);
				break;
			case 'H':
				sequence.push_back(8);
				break;
			case 'I':
				sequence.push_back(9);
				break;
			case 'J':
				sequence.push_back(27);
				break;
			case 'K':
				sequence.push_back(10);
				break;
			case 'L':
				sequence.push_back(11);
				break;
			case 'M':
				sequence.push_back(12);
				break;
			case 'N':
				sequence.push_back(13);
				break;
			case 'O':
				sequence.push_back(26);
				break;
			case 'P':
				sequence.push_back(14);
				break;
			case 'Q':
				sequence.push_back(15);
				break;
			case 'R':
				sequence.push_back(16);
				break;
			case 'S':
				sequence.push_back(17);
				break;
			case 'T':
				sequence.push_back(18);
				break;
			case 'U':
				sequence.push_back(24);
				break;
			case 'V':
				sequence.push_back(19);
				break;
			case 'W':
				sequence.push_back(20);
				break;
			case 'X':
				sequence.push_back(21);
				break;
			case 'Y':
				sequence.push_back(22);
				break;
			case 'Z':
				sequence.push_back(23);
				break;
			case '*':
				sequence.push_back(25);
				break;
			}
		}
	}
}
