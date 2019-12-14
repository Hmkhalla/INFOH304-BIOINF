#define NOMINMAX
#include "Database.h"
#include <algorithm>
using namespace std;

Database::Database(string path)
{

    conversion = {
        {'-', 0}, {'A', 1}, {'B', 2}, {'C', 3}, {'D', 4}, {'E', 5}, {'F', 6}, {'G', 7}, {'H', 8}, {'I', 9}, {'J', 27}, {'K', 10}, {'L', 11}, {'M', 12}, {'N', 13}, {'O', 26}, {'P', 14}, {'Q', 15}, {'R', 16}, {'S', 17}, {'T', 18}, {'U', 24}, {'V', 19}, {'W', 20}, {'X', 21}, {'Y', 22}, {'Z', 23}, {'*', 25} //etc
    };

    readPin(path);
	ifstream psq_stream = ifstream(path + ".psq", ios::binary);
	ifstream phr_stream = ifstream(path + ".phr", ios::binary);
	getBlosumMatrix("BLOSUM62");
	if (psq_stream.is_open() && phr_stream.is_open())
    {
		psq = new uint8_t[Index_seq_table[nb_seq]];
		phr = new uint8_t[Index_head_table[nb_seq]];
		psq_stream.read((char *)&psq[0], sizeof(uint8_t) * (Index_seq_table[nb_seq]));
		phr_stream.read((char *)&phr[0], sizeof(uint8_t) * (Index_head_table[nb_seq]));
		psq_stream.close();
		phr_stream.close();
	}else{
		cerr << "Impossible d'ouvrir le fichier psq ou phr" << endl;
		exit(EXIT_FAILURE);
	}
}

Database::~Database()
{
    delete title;
    delete time;
	delete[] phr;
    delete[] psq;
    delete[] Index_seq_table;
    delete[] Index_head_table;
}

void Database::readPin(string path){
	ifstream pin = ifstream(path + ".pin", ios::binary);
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
    else{
		cerr << "Impossible d'ouvrir le fichier pin" << endl;
		exit(EXIT_FAILURE);
   }
}

int16_t Database::Max(int16_t n1, int16_t n2, int16_t n3, int16_t n4){
	int16_t res = 0;
	if(res < n1) res = n1;
	if(res < n2) res = n2;
	if(res < n3) res = n3;
	if(res < n4) res = n4;
	return res;
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

uint8_t *Database::find_seq(int index, int &nb)
{

    nb = ((int)getIndexSeq(index + 1) - (int)getIndexSeq(index)) - 1;
    return psq + (int)getIndexSeq(index);
}

/*void Database::find_header(string &res, int index)
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
}*/

void Database::exactMatch(char *queryPath)
{/*
    uint8_t *query_seq, temp_seq;
    int nb1, nb2;
    ifstream queryFlux(queryPath, ios::binary);
    query_seq = convertToValue(queryFlux, nb1);
    queryFlux.close();
    for (int i = 0; i < nb_seq; i++)
    {
        //if (getIndexSeq(i + 1) - getIndexSeq(i) - 1 == query_seq.size())
        {
            temp_seq= find_seq(i, nb2);
            for (int i = 0; i <nb1 && ; i ++)
            if (isMatch(query_seq, temp_seq))
            {
                //string header;
                //find_header(header, i);
                cout << "answer here index : " << i << endl;
                //cout << "answer here: " << header << endl;
            }
        }
        temp_seq.clear();
	}
	* */
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
        for (int i = 3; i < 73; i += 3)
        {
			if (line.at(i - 1)== ' ' or line.at(i - 1) == '-') blosumMatrix[conversionCharChar[(int)i / 3]][conversion[line.at(0)]] = ((int)line.at(i) - 48) * ((line.at(i - 1) == '-') ? -1 : 1);
			else blosumMatrix[conversionCharChar[(int)i / 3]][conversion[line.at(0)]] = 11;

        }
    }
    blosumFile.close();
}

int Database::scoring(uint8_t *query, uint32_t queryLen, uint8_t *tmp, uint32_t tmpLen)
{

    int8_t R = 1;
    int8_t open_gp = 11;
    int8_t Q = R+open_gp;
    int8_t blos;
    int16_t F = 0;
    int16_t H_top = 0;
    int16_t H=0;
    int16_t *H_left = new int16_t[queryLen + 1];
    int16_t *E_left = new int16_t[queryLen];
    for (int i = 0; i <=queryLen; i++){
		H_left[i] = 0;
		E_left[i]= 0;
	}

    int16_t max_scoring = 0;
    for (unsigned int j = 1; j <= tmpLen; j++)
    {
        for (unsigned int i = 1; i <= queryLen; i++)
        {
			E_left[i-1] = std::max (H_left[i]-Q, E_left[i-1]-R);
			F = std::max (H_top-Q, F-R);
			blos = blosumMatrix[query[i-1]][tmp[j-1]];
			H = Max(H_left[i-1]+blos, E_left[i-1], F, 0);
			H_left[i-1] = H_top; H_top = H;
			//cout << "(" << i-1 << ", " << j-1 << ") " << "(" << (int )query[i-1] << ", " << (int )tmp[j-1] << ") " << " E : " << E_left[i] << " F : " << F << " Hcurrent : " << H << " Blosum : " << (int) blos <<endl;
			if (H > max_scoring) max_scoring = H;
        }
        F=0;H_top=0;H_left[queryLen] = H;
    }
    delete[] H_left;
    return max_scoring;
}

void Database::notExactMatch(char *queryPath)
{
    uint8_t *query_seq, *temp_seq;
    vector<int> scores;
    int score;
    ifstream queryFlux(queryPath, ios::binary);
    int nb1, nb2;
    query_seq = convertToValue(queryFlux, nb1);
    queryFlux.close();
    scores.resize(nb_seq);
    int mx = 0;
    int maxIndex = 0;
    for (int i = 0; i < nb_seq; i++)
    //for (int i = 2958; i < 2959; i ++)

    {
        if (i % 1000 == 18)
        {
            printf("On est à %d/%d \n", i, nb_seq);
        }
        //printf("On est à %d/%d \n", i, nb_seq);
        temp_seq = find_seq(i, nb2);
        score = scoring(query_seq, nb1, temp_seq, nb2);
        if (score > mx)// && i != 2958)
        {
            mx = score;
            maxIndex = i;
        }
        scores.push_back(score); //min: 32 2968
    }
    sort(scores.begin(), scores.end(), greater <>());
    int j = 0;
    for (const auto &i: scores){
		if (j == 10) break;
		cout << "Score : " << (int) i <<endl;
		j++;
	}/*
    for (int i = 0; i < 10; i++)
		cout << "Score : " << (int) scores[i] << endl;
		*/
    cout << "Max index :" << maxIndex << " Max score : " << mx << endl;
    /*
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
		string header;
		find_header(header, (int)curMaxIndex);
		printf("%s, score %i\n", header, (int)curMax);
	}
	// ===> On créer une liste avec les scores de protéines et et leur index dans les offset*/
}

uint8_t *Database::convertToValue(ifstream &queryFlux, int &nb)
{
    string buffer;
	char residue;
    if (queryFlux.is_open())
    {
        string header;
        getline(queryFlux, header);
        while (queryFlux.get(residue))
        {
            if (residue != 10 && residue != 13)
            {
				
                buffer += residue;
            }
        }
        cout << endl;
        nb = buffer.length();
        uint8_t *sequence = new uint8_t[nb];
        for (int i = 0; i < buffer.length(); i++)
        {
            sequence[i] = conversion[buffer[i]];
        }
        return sequence;
    }
}
