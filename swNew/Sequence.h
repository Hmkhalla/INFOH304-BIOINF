#ifndef DEF_SEQUENCE
#define DEF_SEQUENCE

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

using namespace std;



struct Sequence
{
	uint index;
	int score;
	string name;
	int len;
	uint8_t *sequence;
	
	void convertToValue(string &fasta);

//public:
	Sequence();
	Sequence(string &fasta);
	Sequence(uint8_t *sequence, int len, int index);
	~Sequence();
	const Sequence & operator= ( const Sequence & SequenceOriginal);
	uint8_t operator[] (unsigned int i) const;
	bool operator== (const Sequence &right) const;
	void free_sequence();
	const string getName() const;
	const uint8_t *getSequence() const;
	void setName(string &name);
	void setScore(int nb);
	int getLen() const;
} ;
#endif
