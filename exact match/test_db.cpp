#include "Database.h"

int main (int argc, char *argv[]){
	
	Database * db = new Database(argv[1]);
	db->printDbDescription();
	db->exactMatch(argv[2]);
	delete db;
	return 0; 
}
