#include "Database.h"
using namespace std;

Database::Database(string &path)
{
    readPin(path);
    loadPsq(path);
    loadPhr(path);
}

Database::~Database()
{
	delete[] phr;
    delete[] psq;
    delete[] Index_seq_table;
    delete[] Index_head_table;
}

void Database::loadPhr(string &path)
{
	ifstream phr_stream = ifstream(path + ".phr", ios::binary);
	if (phr_stream.is_open())
    {
		phr = new uint8_t[Index_head_table[nb_seq]];
		phr_stream.read((char *)&phr[0], sizeof(uint8_t) * (Index_head_table[nb_seq]));
		phr_stream.close();
	}else{
		cerr << "Impossible d'ouvrir le fichier phr" << endl;
		exit(EXIT_FAILURE);
	}
}

void Database::loadPsq(string &path)
{
	ifstream psq_stream = ifstream(path + ".psq", ios::binary);
	if (psq_stream.is_open())
    {
		psq = new uint8_t[Index_seq_table[nb_seq]];
		psq_stream.read((char *)&psq[0], sizeof(uint8_t) * (Index_seq_table[nb_seq]));
		psq_stream.close();
	}else{
		cerr << "Impossible d'ouvrir le fichier psq" << endl;
		exit(EXIT_FAILURE);
	}
}


void Database::readPin(string & path)
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
		cerr << "Impossible d'ouvrir le fichier pin" << endl;
		exit(EXIT_FAILURE);
   }
}

void Database::printDbDescription() const
{/*
    cout << "Version :" << version << endl;
    cout << "db_type :" << db_type << endl;
    cout << "title :" << title << endl;
    cout << "Timestamp :" << time << endl;
    cout << "Number of sequences :" << nb_seq << endl;
    cout << "Residue count :" << res_count << endl;
    cout << "Maximum sequence :" << max_seq << endl;
    for (int i = 0; i < 5; i++)
    {
        cout << "Header offset table " << i << " : " << Index_head_table[i] << endl;
    }
    cout << "Header offset table " << nb_seq << " : " << Index_head_table[nb_seq] << endl;
    cout << "Sequence offset table " << nb_seq << " : " << Index_seq_table[nb_seq] << endl;
    for (int i = 0; i < 5; i++)
    {
        cout << "Sequence offset table " << i << " : " << Index_seq_table[i] << endl;
    }*/
    cout << description; 
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

uint8_t *Database::find_seq(int index, int &len) const
{
	len = ((int)getIndexSeq(index + 1) - (int)getIndexSeq(index)) - 1;
	return psq + getIndexSeq(index);
}

void Database::find_header(string &res, int index) const
{
	string title = "";
	string type = "";
	string db = "";
	string id = "";
	int start, end;
	start = getIndexHead(index);
	end = getIndexHead(index + 1);
	bool readingInt = false;
	bool readingLen = false;
	bool readingStr = false;
	uint8_t byte;
	for (int i = 0; i < end; i++){
		byte = phr[start+i];
		if (!readingInt && !readingLen && !readingStr)
		{
			switch (byte)
			{
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
		else if (readingLen)
		{
			unsigned long len = 0;
			if (readingStr)
			{
				if (byte & (1 << 7))
				{
					byte = (byte & ~(1 << 7)) | ((0 << 7) & (1 << 7));
					for (int j = 0; j < byte; j++)
					{
						//phr.read((char *)&buffer, sizeof(uint8_t));
						len = len << 8 | (unsigned long)phr[start+i+j];
					}
					i+=byte;
					//phr.read((char *)&(len), sizeof(byte) * (int)byte);
				}
				else
					len = (int)byte;
				char buffer[len + 1];
				for (int j = 0; j<len; j++){
					buffer[j]=phr[start+i+j];
					//phr.read((char *)&(buffer[0]), sizeof(byte) * len);
				}
				i+=len;
				buffer[len] = 0;
				readingStr = false;
				if (title.empty())
					title = buffer;
				else
					db = buffer;
			}
			else if (readingInt)
			{
				len = (int)byte;
				unsigned long intFound = 0;
				//uint8_t buffer = 0;
				for (int j = 0; j < len; j++)
				{
					//phr.read((char *)&buffer, sizeof(uint8_t));
					intFound = intFound << 8 | (unsigned long)phr[start+i+j];
				}
				i+=len;
				readingInt = false;
				if (id.empty())
					id = to_string((int)intFound);
				else
					type = to_string((int)intFound);
			}
			readingLen = false;
		}
	}
	res = type + "|" + db + "|" + id + " " + title;
}
