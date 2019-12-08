#define NOMINMAX
#include "Database.h"
using namespace std;

Database::Database(string path)
{

	conversion = {
		{'-', 0}, {'A', 1}, {'B', 2}, {'C', 3}, {'D', 4}, {'E', 5}, {'F', 6}, {'G', 7}, {'H', 8}, {'I', 9}, {'J', 27}, {'K', 10}, {'L', 11}, {'M', 12}, {'N', 13}, {'O', 26}, {'P', 14}, {'Q', 15}, {'R', 16}, {'S', 17}, {'T', 18}, {'U', 24}, {'V', 19}, {'W', 20}, {'X', 21}, {'Y', 22}, {'Z', 23}, {'*', 25} //etc
	};

	pin = ifstream(path + ".pin", ios::binary);
	phr = ifstream(path + ".phr", ios::binary);
	psq = ifstream(path + ".psq", ios::binary);
	getBlosumMatrix("BLOSUM62");
	if (pin.is_open())
	{

		uint32_t length;

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

		pin.close();
	}
	else
		cerr << "Impossible d'ouvrir le fichier pin" << endl;
}

Database::~Database()
{
	pin.close();
	phr.close();
	psq.close();
	delete title;
	delete time;
	delete[] Index_seq_table;
	delete[] Index_head_table;
}

void Database::printDbDescription() const
{
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

void Database::find_seq(vector<uint8_t> &bytes, int index)
{

	int len = ((int)getIndexSeq(index + 1) - (int)getIndexSeq(index)) - 1;
	bytes.resize(len);
	if (psq.is_open())
	{
		psq.seekg((int)getIndexSeq(index) * sizeof(uint8_t), ios::beg);
		psq.read((char *)&(bytes[0]), sizeof(uint8_t) * len);
	}
	else
	{
		cerr << "Impossible d'ouvrir le fichier psq" << endl;
	}
}

void Database::find_header(string &res, int index)
{
	if (phr.is_open())
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
		phr.seekg(start * sizeof(byte));
		while (phr.tellg() < end && phr.read((char *)(&byte), sizeof(byte)))
		{
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
				int len = 0;
				if (readingStr)
				{
					if (byte & (1 << 7))
					{
						byte = (byte & ~(1 << 7)) | ((0 << 7) & (1 << 7));
						phr.read((char *)&(len), sizeof(byte) * (int)byte);
					}
					else
						len = (int)byte;
					char buffer[len + 1];
					phr.read((char *)&(buffer[0]), sizeof(byte) * len);
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
					uint8_t buffer = 0;
					for (int j = 0; j < len; j++)
					{
						phr.read((char *)&buffer, sizeof(uint8_t));
						intFound = intFound << 8 | (unsigned long)buffer;
					}
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
	else
		cerr << "Impossible d'ouvrir le fichier" << endl;
}

void Database::exactMatch(char *queryPath)
{
	vector<uint8_t> query_seq, temp_seq;
	ifstream queryFlux(queryPath, ios::binary);
	convertToValue(query_seq, queryFlux);
	queryFlux.close();

	for (int i = 0; i < nb_seq; i++)
	{
		if (getIndexSeq(i + 1) - getIndexSeq(i) - 1 == query_seq.size())
		{
			find_seq(temp_seq, i);
			if (query_seq == temp_seq)
			{
				string header;
				find_header(header, i);
				cout << "index : " << i << endl;
				cout << "answer here: " << header << endl;
			}
		}
		temp_seq.clear();
	}
}

void Database::getBlosumMatrix(string pathBlosum)
{
	ifstream blosumFile;
	blosumFile.open(pathBlosum);
	if (!blosumFile)
	{
		cout << "Unable to open file" << endl;
		exit(1);
	}
	string line;
	map<int, uint8_t> conversionCharChar;
	bool salut = true;
	while (getline(blosumFile, line))
	{
		if (line.at(0) == '#')
			continue;
		if (line.at(0) == ' ')
		{
			for (int i = 3; i < 73; i += 3)
			{
				conversionCharChar.insert({(int)i / 3, conversion[line.at(i)]});
			}
			continue;
		}
		cout << line << endl;
		for (int i = 3; i < 73; i += 3)
		{
			blosumMatrix[conversionCharChar[(int)i / 3]][conversion[line.at(0)]] = ((int)line.at(i) - 48) * ((line.at(i - 1) == '-') ? -1 : 1);
		}
	}

	for (uint8_t i = 0; i < conversion.size(); i++)
	{
		for (uint8_t j = 0; j < conversion.size(); j++)
		{
			if (blosumMatrix[i][j] >= 0)
				cout << " ";
			cout << blosumMatrix[i][j] << " ";
		}
		cout << endl;
	}
	blosumFile.close();
}

uint16_t Database::max(uint16_t n1, uint16_t n2)
{
	return (n1 > n2) ? n1 : n2;
}

int Database::scoring(vector<uint8_t> *seq1, uint32_t nb1, vector<uint8_t> *seq2, uint32_t nb2)
{
	/*
	Here the curV[j][0] = Hij and v2[j] = Hi-1,j
	Here the curV[j][1] = Fij
	We need only two row to have all the information we need
	*/
	uint16_t ext_gap = 1;
	uint16_t open_gp = 11;
	uint16_t sum = 12;
	//vector<int8_t> *v1 = new vector<int8_t>(nb1 + 1, 0);
	uint16_t *curV = new uint16_t[nb2 + 1];
	uint16_t *FijVector = new uint16_t[nb2 + 1];
	for (int i = 0; i < nb2 + 1; i++)
	{
		FijVector[i] = 0;
		curV[i] = 0;
	}

	uint16_t Fij;
	uint16_t Eij;
	uint16_t Hij;
	int16_t max_scoring = 0;
	uint16_t value = 0;
	//printf("le value de depart : %u \n", value);
	int16_t scoreAB;
	//printf( "nb1 = %d et nb2 = %d  et v2 size %d et size seq1 : %d \n", nb1, nb2, v2->size(), seq1->size());

	for (unsigned int i = 1; i <= nb1; i++)
	{
		for (unsigned int j = 1; j <= nb2; j++)
		{
			if (value > Eij + open_gp)
			{
				Eij = value;
				Eij -= sum;
			}
			else if (Eij > ext_gap)
			{
				Eij -= ext_gap;
			}
			else
				Eij = 0;

			//Eij = (Eij > ext_gap) ? Eij - ext_gap : 0;
			//Eij = ((value > Eij + sum) ? value - sum : Eij);
			if (curV[j] > Fij + open_gp)
			{
				Fij = curV[j];
				Fij -= sum;
			}
			else if (Fij > ext_gap)
			{
				Fij -= ext_gap;
			}
			else
				Fij = 0;

			//Fij = (FijVector[j] > ext_gap) ? FijVector[j] - ext_gap : 0;
			//Fij = ((curV[j] > Fij + sum) ? curV[j] - sum : Fij);
			FijVector[j] = Fij;
			scoreAB = blosumMatrix[seq1->at(i - 1)][seq2->at(j - 1)];

			Hij = (curV[j - 1] > -scoreAB) ? curV[j - 1] + scoreAB : 0;

			curV[j - 1] = value;
			value = max(max(Hij, Fij), Eij);
			//tempHij = curV[j];
			FijVector[j] = Fij;

			if (value > max_scoring)
			{
				max_scoring = value;
			}
			//printf("%6u", value);
		}
		//*curV = vector<uint16_t>(nb2 + 1, 0);
		//printf("\n");
		value = 0;
		Eij = 0;
	}
	//printf("max scoring = %d", max_scoring);
	delete[] curV;
	delete[] FijVector;
	return max_scoring;
}

void Database::notExactMatch(char *queryPath)
{
	vector<uint8_t> query_seq, temp_seq;
	vector<uint8_t> scores;
	int score;
	ifstream queryFlux(queryPath, ios::binary);
	convertToValue(query_seq, queryFlux);
	queryFlux.close();
	scores.resize(nb_seq);
	int mx = 0;
	int maxIndex = 0;
	for (int i = 0; i < nb_seq; i++)
	{
		if (i % 1000 == 18)
		{
			printf("On est à %d/%d \n", i, nb_seq);
		}
		//printf("On est à %d/%d \n", i, nb_seq);
		find_seq(temp_seq, i);

		score = scoring(&query_seq, query_seq.size(), &temp_seq, temp_seq.size());
		if (score > mx)
		{
			mx = score;
			maxIndex = i;
		}
		//scores.push_back(score); //min: 32
	}
	cout << "Max index :" << maxIndex << endl;

	// ===> On créer une liste avec les scores de protéines et et leur index dans les offset*/
}

void Database::convertToValue(vector<uint8_t> &sequence, ifstream &queryFlux)
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
			sequence.push_back(conversion[buffer[i]]);
			/*
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
			}*/
		}
	}
}
