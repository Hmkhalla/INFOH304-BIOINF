#include<iostream>
#include<fstream>
using namespace std;


int main(int argc, char *argv[])
{
	ifstream in (argv[1], ios::binary);
	//ofstream out ("out.txt");
	if( in.is_open() )
	{	
		int nbStrRead = 0;
		int nbIntRead = 0;
		int length=0;
		bool readingInt =false;
		bool readingLen =false;
		bool readingStr =false;
		uint8_t byte;
		while( in.read((char *) (&byte), sizeof(byte) )){
			cout << std::hex << (int) byte << endl;
			if ( ! readingInt && ! readingLen && !readingStr){
				switch(byte) {
				  case 0x1a:
					// String
					readingStr = true;
					readingLen = true;
					cout << "str detected" << endl;
					break;
				  case 0x02:
					// int
					readingInt = true;
					readingLen = true;
					cout << "int detected" << endl;
					break;
				}
			}
			else if (readingLen){
				
				if (readingStr){
					if (byte & (1<<7)){
						cout << "Long lentgh - Str" << endl;
						readingStr = false;
						readingLen = false;
					}else{
						char title[(int) byte +1];
						for (int i=0; i < (int) byte; i++){
							in.read((char *) &(title[i]), sizeof(byte));
						}
						title[(int) byte] = 0;
						readingStr = false;
						readingLen = false;
						cout << title << endl;
						cout << "end str" <<endl;
					}
				}else if (readingInt){
					cout << "read int" << endl;
					readingInt = false;
					readingLen = false;
				}
			}
			
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
