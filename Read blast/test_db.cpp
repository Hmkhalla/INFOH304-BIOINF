#include "Database.h"

int main (int argc, char *argv[]){
	cout<<"1" << endl;
	Database * db = new Database(argv[1]);
	cout<<"2" << endl;
	//const uint32_t *a = db->getSeqTable();
	db->printDbDescription();
	cout << db->prot_table[15678].getProteinDescription() << endl;
	//db->prot_table[1].printVectorSequence();
	//for (const uint8_t *i = seq; (int)*i!=0; i++) cout << (int) *i <<endl;
	delete db;
	return 0; 
}
