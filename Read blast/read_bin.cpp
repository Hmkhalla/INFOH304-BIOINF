#include<iostream>
#include<fstream>
using namespace std;


int main(int argc, char *argv[])
{
	ifstream in ("P00533.fasta.psq", ios::binary);
	ofstream out ("out.txt");
	if( in.is_open() )
	{
		int8_t x;
		while( in.read((char *) (&x), sizeof(x) )){
			cout << (int) x << endl;
			out << (int) x << endl;
		}
		in.close();
		out.close();
		//cout << "longueur " << l.getNombreNoeuds() << endl;
		//l.imprimeListe(out);
	}
	else
		cout << "Impossible d'ouvrir le fichier" << endl;
	return 0;
}
