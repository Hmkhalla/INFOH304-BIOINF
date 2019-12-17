#include "Algorithm.h"

using namespace std;

Algorithm::Algorithm(string dbPath, string fasta, string pathBlosum, int core, int extension_gap, int open_gap)
{
	conversion = {
        {'-', 0}, {'A', 1}, {'B', 2}, {'C', 3}, {'D', 4}, {'E', 5}, {'F', 6}, {'G', 7}, {'H', 8}, {'I', 9}, {'J', 27}, {'K', 10}, {'L', 11}, {'M', 12}, {'N', 13}, {'O', 26}, {'P', 14}, {'Q', 15}, {'R', 16}, {'S', 17}, {'T', 18}, {'U', 24}, {'V', 19}, {'W', 20}, {'X', 21}, {'Y', 22}, {'Z', 23}, {'*', 25} //etc
    };
	R=extension_gap;
	Q=extension_gap + open_gap;
	db = new Database(dbPath);
	nb_seq = db->getNbSeq();
	query = Sequence(fasta);
	initBlosumMatrix(pathBlosum);
	nb_thread= core; 
	threads = new thread[nb_thread];
	seqArray = new Sequence[nb_seq];
	db->printDbDescription();
}

Algorithm::~Algorithm()
{
	delete db;
	delete []threads;
	delete []seqArray;
	query.free_sequence();
}
	
int Algorithm::Max(int n1, int n2, int n3) const{
	int res = 0;
	if(res < n1) res = n1;
	if(res < n2) res = n2;
	if(res < n3) res = n3;
	return res;
}

void Algorithm::initBlosumMatrix(string &pathBlosum)
{
    ifstream blosumFile;
    blosumFile.open(pathBlosum);
    if (!blosumFile)
    {
        cout << "Unable to open file" << endl;
        exit(1);
    }
    for (int i = 0; i<28; i++) {
		for (int j = 0; j<28; j++) {
			blosumMatrix[i][j] = 0; 
		}
	}
    string line;
    map<int, uint8_t> conversionChar;
    while (getline(blosumFile, line))
    {
        if (line.at(0) == '#')
            continue;
        if (line.at(0) == ' ')
        {
            for (int i = 3; i < 73; i += 3)
            {
                conversionChar.insert({(int)i / 3, conversion[line.at(i)]});
            }
            continue;
        }
        for (int i = 3; i < 73; i += 3)
        {
			if (line.at(i - 1)== ' ' or line.at(i - 1) == '-') blosumMatrix[conversionChar[(int)i / 3]][conversion[line.at(0)]] = ((int)line.at(i) - 48) * ((line.at(i - 1) == '-') ? -1 : 1);
			else blosumMatrix[conversionChar[(int)i / 3]][conversion[line.at(0)]] = 11;

        }
    }
    blosumFile.close();
}

int Algorithm::scoring(Sequence& tmp) const
{
    int blos;
    int F = 0;
    int H_top = 0;
    int H=0;
    int queryLen = query.getLen(); int tmpLen = tmp.getLen();
    int *H_left = new int[queryLen + 1];
    int *E_left = new int[queryLen + 1];
    for (int i = 0; i <=queryLen; i++){
		H_left[i] = 0;
		E_left[i]= 0;
	}
    int max_scoring = 0;
    for (unsigned int j = 1; j <= tmpLen; j++)
    {
        for (unsigned int i = 1; i <= queryLen; i++)
        {
			E_left[i] = std::max(H_left[i]-Q, E_left[i]-R);
			F = std::max (H_top-Q, F-R);
			blos = blosumMatrix[query.sequence[i-1]][tmp.sequence[j-1]];
			H = Max(H_left[i-1]+blos, E_left[i], F);
			H_left[i-1] = H_top; H_top = H;
			//if (i == j)
			//cout << "(" << i-1 << ", " << j-1 << ") " << "(" << (int )query[i-1] << ", " << (int )tmp[j-1] << ") " << " E : " << E_left[i] << " F : " << F << " Hcurrent : " << H << " Blosum : " << (int) blos <<endl;
			if (H > max_scoring) max_scoring = H;
        }
        F=0;H_top=0;H_left[queryLen] = H;
    }
    delete[] H_left;
    delete[] E_left;
    return max_scoring;
}

void Algorithm::startMultithread(){
	for(int i = 0; i<nb_thread; i++){
	    threads[i]=std::thread(&Algorithm::swAlgo, ref(*this), i); 
	}
	for (int i = 0; i<nb_thread; i++) threads[i].join();
}

void Algorithm::showResult(){
	std::sort(seqArray, seqArray+nb_seq, [](Sequence &a, Sequence &b) {
			return a.score > b.score;
		});
    for (int i = 0; i<10;i++){
		string header;
		db->find_header(header, seqArray[i].index);
		cout << "Score : " << (int) seqArray[i].score <<endl << header<< endl <<endl;
	}
}

void Algorithm::swAlgo(int start)
{
    Sequence temp_seq;
    int score;
    int tmp_len;
    int a;
    uint8_t* point;
    for (int i = start; i < nb_seq; i+=nb_thread)
    {
        /*if (i % 1000 == 18)
        {
            printf("On est à %d/%d \n", i, nb_seq);
        }*/
        point=db->getSeq(i, a);
        temp_seq = Sequence(point, a, i);
        score = scoring(temp_seq);
        temp_seq.setScore(score);
        seqArray[i]=temp_seq; //min: 32 2968
    }
}

void Algorithm::exactMatch()
{
    Sequence temp_seq;
    int tmp_len=0;
    for (int i = 0; i < nb_seq; i++)
    {
        if (db->getLenSeq(i) == query.getLen())
        {
			uint8_t* seq = db->getSeq(i, tmp_len);
			temp_seq= Sequence(seq, tmp_len, i);
			cout <<"admissible "<<endl;
            if (query == temp_seq)
            {
                string header;
                db->find_header(header, i);
                cout << "index : " << i << endl;
                cout << "answer here: " << header << endl;
            }
        }
    }
}


