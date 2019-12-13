#include "Protein.h"
using namespace std;

uint8_t *Protein::convertToValue(string m_proteinSequence) //why int &nb?
{
    string buffer;
    buffer = m_proteinSequence;
    int nb = buffer.length();
    uint8_t *sequence = new uint8_t[nb];
    for (int i = 0; i < buffer.length(); i++)
    {
        sequence[i] = m_conversion[buffer[i]];
    }
    
        return sequence;
}

string Protein::getProteinDescription()
{
	return m_proteinDescription;
}

string Protein::getProteinSequence()
{
	return m_proteinSequence;
}

map<char, uint8_t> Protein::getConversion()
{
	return m_conversion;
}
	


uint8_t *Protein::getSequenceConverted()
{
	return m_sequenceConverted;
}

Protein::Protein(string fasta)
{
	ifstream queryFlux(fasta);
	if (queryFlux.is_open())
	{
		string line;
		while (getline(queryFlux, line).good())
		{
			if (line[0] == '>')
			{
				m_proteinDescription = line.substr(1);
			}
			else
			{
				m_proteinSequence += line;
			}
		}
		m_conversion = 
{
  {'-', 0}, {'A', 1}, {'B', 2}, {'C', 3}, {'D', 4}, {'E', 5}, {'F', 6}, {'G', 7}, {'H', 8}, {'I', 9}, {'J', 27}, {'K', 10}, {'L', 11}, {'M', 12}, {'N', 13}, {'O', 26}, {'P', 14}, {'Q', 15}, {'R', 16}, {'S', 17}, {'T', 18}, {'U', 24}, {'V', 19}, {'W', 20}, {'X', 21}, {'Y', 22}, {'Z', 23}, {'*', 25} //etc
};
		m_sequenceConverted=convertToValue(m_proteinSequence);
	}
	else
	{
		cout << "Error: Impossible to open the file." << endl;
	}
}

Protein::~Protein()
{
	delete m_sequenceConverted;
}





	




