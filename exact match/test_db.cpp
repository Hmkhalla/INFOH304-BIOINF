#include "Database.h"
using namespace std;

int main (int argc, char *argv[]){
	
	for(int i=1;i<=2;i++){
	    string file= argv[i];
	    if (file.substr(file.find_last_of(".") + 1) != "fasta"){
		    std::cerr << "Usage: " << argv[i] <<"Enter: fasta file" << std::endl;
		    return 1; }
	}
	
	Database * db = new Database(argv[1]);
	db->printDbDescription();
	db->exactMatch(argv[2]);
	delete db;
	return 0; 
}
