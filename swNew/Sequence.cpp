#include "Sequence.h"
using namespace std;

Sequence::Sequence() 
{
	index=0;
	score=0;
	len=0;
	sequence=NULL;
}

Sequence::Sequence(string &fasta)  //constructor for the query sequence
{	
	convertToValue(fasta);
	cout << "Query file name: " << fasta << endl << "Query length: " << len << " residues" << endl << "Query description: " << name << endl << endl;
}

Sequence::Sequence(uint8_t *table, int ln, int i) //constructor for the template sequence from the database
{
	sequence=table;
	len=ln;
	index=i;
}

const Sequence & Sequence::operator= ( const Sequence & SequenceOriginal) //redefines the assignment operator for sequences
{
	if ( this != &SequenceOriginal ){
		index=SequenceOriginal.index;
		score=SequenceOriginal.score;
		name=SequenceOriginal.name;
		len=SequenceOriginal.len;
		sequence=SequenceOriginal.sequence;
		return *this;
	}
}

bool Sequence::operator== (const Sequence &right) const //redefines the equality operator for comparing sequences
{
	int i;
	if (len != right.len) return false;
	for (i = 0; sequence[i]==right.sequence[i] && i<len; i++);
	if (i != len){
		return false;
	 }
	return true;
}

	
void Sequence::free_sequence()
{
	delete [] sequence;
}

Sequence::~Sequence()
{
}

void Sequence::convertToValue(string &fasta) //converts each residue of a protein sequence,encoded as a char, to its value in int
{
	map<char, uint8_t> conversion = {
        {'-', 0}, {'A', 1}, {'B', 2}, {'C', 3}, {'D', 4}, {'E', 5}, {'F', 6}, {'G', 7}, {'H', 8}, {'I', 9}, {'J', 27}, {'K', 10}, {'L', 11}, {'M', 12}, {'N', 13}, {'O', 26}, {'P', 14}, {'Q', 15}, {'R', 16}, {'S', 17}, {'T', 18}, {'U', 24}, {'V', 19}, {'W', 20}, {'X', 21}, {'Y', 22}, {'Z', 23}, {'*', 25} //etc
    };
    string buffer;
	char residue;
	ifstream queryFlux(fasta, ios::binary);
    if (queryFlux.is_open())
    {
        getline(queryFlux, name);
        while (queryFlux.get(residue))
        {
            if (residue != 10 && residue != 13)
            {
                buffer += residue;
            }
        }
        len = buffer.length();
        sequence = new uint8_t[len];
        for (int i = 0; i < len; i++)
        {
            sequence[i] = conversion[buffer[i]];
        }
    }else{
		cerr << "Impossible to open the file fasta" << endl;
		exit(EXIT_FAILURE);
	}
}

const string Sequence::getName() const{
	return name;}
	
const uint8_t *Sequence::getSequence() const{
	return sequence;}
	
void Sequence::setName(string &str){
	name=str;}

void Sequence::setScore(int nb){
	score=nb;}
	
int Sequence::getLen() const{
	return len;}
