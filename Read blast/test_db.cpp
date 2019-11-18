#include "Database.h"

int main (int argc, char *argv[]){
	Database * db = new Database(argv[1]);
	
	const uint32_t *a = db->getSeqTable();
	db->printDbDescription();
	const uint8_t * seq = db->find_seq(2);
	for (const uint8_t *i = seq; (int)*i!=0; i++) cout << (int) *i <<endl;
	delete db;
	return 0; 
}
