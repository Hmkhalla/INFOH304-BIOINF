#include<iostream>
#include<fstream>
using namespace std;

int main(int argc, char *argv[])
{
	ifstream in (argv[1], ios::binary);
	//ofstream out ("out_pin.txt");
	if( in.is_open() )
	{
		int32_t Version;
		int32_t db_type;
		int32_t nb_seq;
		int64_t res_c;
		int32_t max_seq;
		uint32_t title_length;
		
		in.read((char *) (&Version), sizeof(Version) );Version=__bswap_32(Version);
		in.read((char *) (&db_type), sizeof(db_type) );db_type=__bswap_32(db_type);
		in.read((char *) (&title_length), sizeof(title_length) );title_length = __bswap_32(title_length);
		char title[title_length];
		in.read((char *) (&title), sizeof(char)*title_length); title[title_length]=0;
		uint32_t time_length;
		in.read((char *) (&time_length), sizeof(time_length) );time_length=__bswap_32(time_length);
		char time[time_length];
		in.read((char *) (&time), sizeof(char)*time_length); time[time_length]=0;
		in.read((char *) (&nb_seq), sizeof(nb_seq) );nb_seq=__bswap_32(nb_seq);
		in.read((char *) (&res_c), sizeof(res_c) );
		in.read((char *) (&max_seq), sizeof(max_seq) );max_seq=__bswap_32(max_seq);
		int32_t header_table_offset[nb_seq+1]; int32_t seq_table_offset[nb_seq+1];
		
		in.read((char *) (&header_table_offset), sizeof(int32_t)*(nb_seq+1) );//header_table_offset[i]=__bswap_32(header_table_offset[i]);
		in.read((char *) (&seq_table_offset), sizeof(int32_t)*(nb_seq+1));//seq_table_offset[i]=__bswap_32(seq_table_offset[i]);
		for (int i = 0; i <nb_seq+1; i++){
			header_table_offset[i]=__bswap_32(header_table_offset[i]);
			seq_table_offset[i]=__bswap_32(seq_table_offset[i]);
		}
		in.close();
		//out.close();
		cout << "Version :" << Version << endl;
		cout << "db_type :" << db_type << endl;
		cout << "title_length :" << title_length << endl;
		cout << "title :" << title << endl;
		cout << "Timestamp length :" << time_length << endl;
		cout << "Timestamp :" << time << endl;
		cout << "Number of sequences :" << nb_seq << endl;
		cout << "Residue count :" << res_c << endl;
		cout << "Maximum sequence :" << max_seq << endl;
		for (int i =0; i< 5; i++){
			cout << "Header offset table " << i << " : "<< header_table_offset[i] << endl;
		}
		int i = nb_seq;
		cout << "Header offset table " << i << " : "<< header_table_offset[i] << endl;
		cout << "Sequence offset table " << i << " : "<< seq_table_offset[i] << endl;
		for (int i =0; i< 5; i++){
			cout << "Sequence offset table " << i << " : "<< seq_table_offset[i] << endl;
		}
	}
	else
		cout << "Impossible d'ouvrir le fichier" << endl;
	return 0;
}
