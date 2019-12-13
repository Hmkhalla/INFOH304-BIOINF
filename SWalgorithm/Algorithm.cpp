#define NOMINMAX 
#include "Algorithm.h"

using namespace std;

Algorithm::Algorithm(string dbPath, string fasta, string pathBlosum)
{
	db = new Database(dbPath);
	db->printDbDescription();
	pr = new Protein(fasta);
	
	map<char, uint8_t> conversion = pr->getConversion();
	
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
            m_blosumMatrix[conversionCharChar[(int)i / 3]][conversion[line.at(0)]] = ((int)line.at(i) - 48) * ((line.at(i - 1) == '-') ? -1 : 1);
        }
    }
    blosumFile.close();
    
    m_extGap = 1;
    m_openGap = 11;
}

Algorithm::~Algorithm()
{
	delete db;
	delete pr;
}
	

int8_t Algorithm::getBlosumMatrix()
{
    return m_blosumMatrix[28][28];
}

int Algorithm::scoring(uint8_t *seq1, uint32_t nb1, uint8_t *seq2, uint32_t nb2)
{
    uint16_t sum = m_extGap + m_openGap;
    uint16_t *curV = new uint16_t[nb2 + 1];
    uint16_t *FijVector = new uint16_t[nb2 + 1];
    for (int i = 0; i < nb2 + 1; i++)
    {
        FijVector[i] = 0;
        curV[i] = 0;
    }

    uint16_t Fij = 0;
    uint16_t Eij = 0;
    int16_t Hij;
    int16_t max_scoring = 0;
    uint16_t value = 0;
    //printf("le value de depart : %u \n", value);
    int16_t scoreAB = 0;
    //printf( "nb1 = %d et nb2 = %d  et v2 size %d et size seq1 : %d \n", nb1, nb2, v2->size(), seq1->size());

    for (unsigned int i = 1; i <= nb1; i++)
    {
        for (unsigned int j = 1; j <= nb2; j++)
        {

            if (value > Eij + m_openGap && value >= sum)
            {
                Eij = value;
                Eij -= sum;
            }
            else if (Eij >= m_extGap)
            {
                Eij -= m_extGap;
            }
            else
                Eij = 0;

            if (curV[j] > FijVector[j] + m_openGap && curV[j] >= sum)
            {
                Fij = curV[j];
                Fij -= sum;
            }
            else if (FijVector[j] >= m_extGap)
            {
                Fij = FijVector[j] - m_extGap;
            }
            else
                Fij = 0;
            /*
            Eij = (Eij > ext_gap) ? Eij - ext_gap : 0;
            Eij = ((value > Eij + sum) ? value - sum : Eij);
            Fij = (FijVector[j] > ext_gap) ? FijVector[j] - ext_gap : 0;
            Fij = ((curV[j] > Fij + sum) ? curV[j] - sum : Fij);
			*/
            FijVector[j] = Fij;
            //scoreAB = blosumMatrix[seq1[i - 1]][seq2[j - 1]];

            Hij = curV[j - 1];
            Hij += m_blosumMatrix[seq1[i - 1]][seq2[j - 1]];

            curV[j - 1] = value;
            //value = max(max(Hij, Eij), Fij);
            if (Hij > Eij)
                value = Hij;
            else
                value = Eij;
            if (Fij > value)
                value = Fij;

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

void Algorithm::notExactMatch()
{
    uint8_t *query_seq, *temp_seq;
    vector<uint8_t> scores;
    int score;
    //ifstream queryFlux(queryPath, ios::binary);
    int nb1, nb2;
    query_seq = pr->getSequenceConverted();
    nb1 = pr->getProteinSequence().length();
    uint32_t nbSequence=db->getNbSeq();
    scores.resize(nbSequence);
    int mx = 0;
    int maxIndex = 0;
    for (int i = 0; i < nbSequence; i++)
    //for (int i = 2958; i < 2959; i ++)

    {
        if (i % 1000 == 18)
        {
            printf("On est à %d/%d \n", i, nbSequence);
        }
        //printf("On est à %d/%d \n", i, nb_seq);
        temp_seq = db->find_seq(i, nb2);

        score = scoring(query_seq, nb1, temp_seq, nb2);
        if (score > mx && i != 2958)
        {
            mx = score;
            maxIndex = i;
        }
        //scores.push_back(score); //min: 32 2968
    }
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


/*void Algorithm::exactMatch(char *queryPath)
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
}*/


