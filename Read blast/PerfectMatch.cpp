#include <iostream>
//#include <fstream>
#include <vector>
#include "Database.h"
#include "Protein.h"
using namespace std;

int main(int argc, char *argv[])
/*
    Function of the perfect match: 
        -> first arg : file .pin
        -> second arg : file .phr
        -> third arg : file .psq
        -> 4th arg : file .fasta
*/
{
    Database *database = new Database(argv[1]);
    Protein prot = Protein(argv[4]); //create a prot class
    const int32_t nb_seq = database->prot_table.size();
    // vector<int32_t> seq_offset = database->get_prot();
    int cur_size;
    int seq_size = prot.m_vectorSequence.size();
    vector<int32_t> ffo; // first filter offsets

    for (int i = 0; i < nb_seq - 1; i++)
    {
        database->prot_table[i].m_vectorSequence.size();
        if (cur_size == seq_size)
        {
            ffo.push_back(i);
        }
    }
    int32_t cur_length;
    char *cur_seq;
    vector<int32_t> lfo;
    for (int i = 0; i < ffo.size(); i++)
    {
        if (database->prot_table[ffo[i]].m_vectorSequence == prot.m_vectorSequence)
        {
            lfo.push_back(ffo[i]);
        } // get_seq_from_offset(FILE file, int32_t offset)
    }
    for (int i = 0; i < lfo.size(); i++)
    {
        cout << database->prot_table[lfo[i]].getProteinDescription() << endl;
    }
}
