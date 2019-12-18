#include "Database.h"
using namespace std;

Database::Database(string &path)
{
    readPin(path);
    loadPsq(path);
    loadPhr(path);
    cout << description<<endl; 
}

Database::~Database()
{
	delete[] phr;
    delete[] psq;
    delete[] Index_seq_table;
    delete[] Index_head_table;
}

void Database::loadPhr(string &path)    //Loading of the Header file into memory 
{
	ifstream phr_stream = ifstream(path + ".phr", ios::binary);
	if (phr_stream.is_open())
    {
		phr = new uint8_t[Index_head_table[nb_seq]];
		phr_stream.read((char *)&phr[0], sizeof(uint8_t) * (Index_head_table[nb_seq]));
		phr_stream.close();
	}else{
		cerr << "Impossible to open the file phr" << endl;
		exit(EXIT_FAILURE);
	}
}

void Database::loadPsq(string &path)   //Loading of the Sequence file into memory
{
	ifstream psq_stream = ifstream(path + ".psq", ios::binary);
	if (psq_stream.is_open())
    {
		psq = new uint8_t[Index_seq_table[nb_seq]];
		psq_stream.read((char *)&psq[0], sizeof(uint8_t) * (Index_seq_table[nb_seq]));
		psq_stream.close();
	}else{
		cerr << "Impossible to open the file psq" << endl;
		exit(EXIT_FAILURE);
	}
}


void Database::readPin(string & path) //Parsing of the Index file and storing the database description
{
	ifstream pin = ifstream(path + ".pin", ios::binary);
	if (pin.is_open())
    {
		std::stringstream stream;
        uint32_t length, version, db_type, max_seq;
		uint64_t res_count;
		char *title;
		char *time;

        pin.read((char *)(&version), sizeof(version));
        version = __bswap_32(version);
        pin.read((char *)(&db_type), sizeof(db_type));
        db_type = __bswap_32(db_type);

        pin.read((char *)(&length), sizeof(length));
        length = __bswap_32(length);
        title = new char[length + 1];
        pin.read((char *)&title[0], sizeof(char) * length);
        title[length] = 0;

        pin.read((char *)&length, sizeof(length));
        length = __bswap_32(length);
        time = new char[length + 1];
        pin.read((char *)&time[0], sizeof(char) * length);
        time[length] = 0;

        pin.read((char *)(&nb_seq), sizeof(nb_seq));
        nb_seq = __bswap_32(nb_seq);
        pin.read((char *)(&res_count), sizeof(res_count));
        pin.read((char *)(&max_seq), sizeof(max_seq));
        max_seq = __bswap_32(max_seq);

        Index_seq_table = new uint32_t[nb_seq + 1];
        Index_head_table = new uint32_t[nb_seq + 1];

        pin.read((char *)&Index_head_table[0], sizeof(int32_t) * (nb_seq + 1));
        pin.read((char *)&Index_seq_table[0], sizeof(int32_t) * (nb_seq + 1));
        for (int i = 0; i <= nb_seq; i++)
        {
            Index_seq_table[i] = __bswap_32(Index_seq_table[i]);
            Index_head_table[i] = __bswap_32(Index_head_table[i]);
        }
        
        stream << "Database file: " << path << endl;
		stream << "Database title: " << path << endl;
        stream << "Database time: " << time << endl;
        stream << "Database size: " << res_count << " in " << nb_seq << " sequences" << endl;
        stream << "Longest db seq: " << max_seq << endl;
		description = stream.str();
        pin.close();
        delete time; delete title;
    }
    else{
		cerr << "Impossible to open the file pin" << endl;
		exit(EXIT_FAILURE);
   }
}

uint32_t Database::getIndexSeq(int index) const
{
    if (index > nb_seq)
    {
        cerr << "Seq index " << index << " out of range" << endl;
        return 0;
    }
    return Index_seq_table[index];
}

uint32_t Database::getIndexHead(int index) const
{
    if (index > nb_seq)
    {
        cerr << "Head index " << index << " out of range" << endl;
        return 0;
    }
    return Index_head_table[index];
}

uint32_t Database::getNbSeq() const
{
	return nb_seq;
}

uint32_t Database::getLenSeq(int index) const{
	return getIndexSeq(index + 1) - getIndexSeq(index) - 1;
}

uint8_t *Database::getSeq(int index, int &len) const
{
	len = ((int)getIndexSeq(index + 1) - (int)getIndexSeq(index)) - 1;
	return psq + getIndexSeq(index);
}

uint8_t *Database::getHeader(int index, int &len) const
{
	len = ((int)getIndexHead(index + 1) - (int)getIndexHead(index)) - 1;
	return phr + getIndexHead(index);
}
	

void Database::find_header(string &description, int index){
	int headerSize;
	uint8_t * header = getHeader(index, headerSize);
	string title, db, type, id;

	int intDone = 0;
	int strDone = 0;
	int len = 0; //length of int or str (in bytes)
	stringstream strBuffer[2];
	stringstream intBuffer[2];
	
	for (int i = 0; i<headerSize; i++){ //reading every bytes of the header
		if(header[i] == 0x1a){//READING STRING
			//reading length 
			i++;
			if(header[i] & (1 << 7)){ //MSB ON
				int a = (header[i] & ~(1 << 7)) | ((0 << 7) & (1 << 7)); //length of the string is stored in "numberOfBytesString" bytes
				for (int j = 0; j < a; j++){ 
					i++;
					len = len << 8 | (int) header[i];
				}
			}
			else{
				//MSB OFF
				len = (int) header[i]; //length of the string is stored int 1 byte
			}
			
			//reading string
			for (int j = 0; j < len; j++){
				i++;
				strBuffer[strDone] << (char) header[i];
			}
			strDone++;
			
		}
		else if(header[i] == 0x02){// READING INT
			//reading length
			i++;
			len = (int) header[i];
 			int intFound= 0;
			//reading int
			for (int j = 0; j < len; j++){ 
				i++;
				intFound = intFound << 8 | (int)header[i];
			}
		    intBuffer[intDone] << intFound;
			
			intDone++;
		}
	}
	description = intBuffer[1].str()+"|"+strBuffer[1].str()+"|"+intBuffer[0].str()+" "+strBuffer[0].str();	
}
