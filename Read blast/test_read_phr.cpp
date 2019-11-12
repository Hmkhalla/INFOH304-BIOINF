#include<iostream>
#include<fstream>
using namespace std;


int main(int argc, char *argv[])
{
	ifstream in (argv[1], ios::binary);
	//ofstream out ("out.txt");
	if( in.is_open() )
	{	
		int index=0;
		int8_t x;
		while( in.read((char *) (&x), sizeof(x) )){
			//cout << (int) x << endl;
			cout << std::hex << (int) x << endl;
		}
		cout << endl;
		in.close();
		//out.close();
		//cout << "longueur " << l.getNombreNoeuds() << endl;
		//l.imprimeListe(out);
	}
	else
		cout << "Impossible d'ouvrir le fichier" << endl;
	return 0;
}
