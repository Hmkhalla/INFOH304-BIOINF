#include "Algorithm.h"
using namespace std;

int main(int argc, char *argv[])
{
	for (int i = 1; i <= 2; i++)
	{
		string file = argv[i];
		if (file.substr(file.find_last_of(".") + 1) != "fasta")
		{
			std::cerr << "Usage: " << argv[i] << "Enter: fasta file" << std::endl;
			return 1;
		}
	}

	Algorithm *al = new Algorithm(argv[1],argv[2],argv[3]);
	//al->exactMatch();
	al->notExactMatch();

	delete al;
	return 0;
}
	


