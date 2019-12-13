#define NOMINMAX
#include "Database.h"
using namespace std;

Database::Database(string dbPath)
{
    pin = ifstream(dbPath + ".pin", ios::binary);
    phr = ifstream(dbPath + ".phr", ios::binary);
    psq = ifstream(dbPath + ".psq", ios::binary);
 
    if (pin.is_open())
    {
        uint32_t length;

        pin.read((char *)(&version), sizeof(version));
        version = __bswap_32(version);
        pin.read((char *)(&db_type), sizeof(db_type));
        db_type = __bswap_32(db_type);

        pin.read((char *)(&length), sizeof(length));
        length = __bswap_32(length);
        title = new char[length + 1];
        pin.read((char *)&title[0], sizeof(char) * length);
        title[length] = 0;

        pin.read((char *)&length, sizeof(length));
        length = __bswap_32(length);
        time = new char[length + 1];
        pin.read((char *)&time[0], sizeof(char) * length);
        time[length] = 0;

        pin.read((char *)(&nb_seq), sizeof(nb_seq));
        nb_seq = __bswap_32(nb_seq);
        pin.read((char *)(&res_count), sizeof(res_count));
        pin.read((char *)(&max_seq), sizeof(max_seq));
        max_seq = __bswap_32(max_seq);

        Index_seq_table = new uint32_t[nb_seq + 1];
        Index_head_table = new uint32_t[nb_seq + 1];

        pin.read((char *)&Index_head_table[0], sizeof(int32_t) * (nb_seq + 1));
        pin.read((char *)&Index_seq_table[0], sizeof(int32_t) * (nb_seq + 1));
        for (int i = 0; i <= nb_seq; i++)
        {
            Index_seq_table[i] = __bswap_32(Index_seq_table[i]);
            Index_head_table[i] = __bswap_32(Index_head_table[i]);
        }

        pin.close();
    }
    else
        cerr << "Impossible d'ouvrir le fichier pin" << endl;
}

Database::~Database()
{
    pin.close();
    phr.close();
    psq.close();
    delete title;
    delete time;
    delete[] Index_seq_table;
    delete[] Index_head_table;
}

void Database::printDbDescription() const
{
    cout << "Version :" << version << endl;
    cout << "db_type :" << db_type << endl;
    cout << "title :" << title << endl;
    cout << "Timestamp :" << time << endl;
    cout << "Number of sequences :" << nb_seq << endl;
    cout << "Residue count :" << res_count << endl;
    cout << "Maximum sequence :" << max_seq << endl;
    for (int i = 0; i < 5; i++)
    {
        cout << "Header offset table " << i << " : " << Index_head_table[i] << endl;
    }
    cout << "Header offset table " << nb_seq << " : " << Index_head_table[nb_seq] << endl;
    cout << "Sequence offset table " << nb_seq << " : " << Index_seq_table[nb_seq] << endl;
    for (int i = 0; i < 5; i++)
    {
        cout << "Sequence offset table " << i << " : " << Index_seq_table[i] << endl;
    }
}

uint32_t Database::getIndexSeq(int index) const
{
    if (index > nb_seq)
    {
        cerr << "Seq index " << index << " out of range" << endl;
        return 0;
    }
    return Index_seq_table[index];
}

uint32_t Database::getIndexHead(int index) const
{
    if (index > nb_seq)
    {
        cerr << "Head index " << index << " out of range" << endl;
        return 0;
    }
    return Index_head_table[index];
}

uint32_t Database::getNbSeq()
{
	return nb_seq;
}

uint8_t *Database::find_seq(int index, int &nb)
{

    nb = ((int)getIndexSeq(index + 1) - (int)getIndexSeq(index)) - 1;
    uint8_t *bytes = new uint8_t[nb];
    if (psq.is_open())
    {
        psq.seekg((int)getIndexSeq(index) * sizeof(uint8_t), ios::beg);
        psq.read((char *)&(bytes[0]), sizeof(uint8_t) * nb);
        return bytes;
    }
    else
    {
        cerr << "Impossible d'ouvrir le fichier psq" << endl;
    }
}

void Database::find_header(string &res, int index)
{
    if (phr.is_open())
    {
        string title = "";
        string type = "";
        string db = "";
        string id = "";
        int start, end;
        start = getIndexHead(index);
        end = getIndexHead(index + 1);
        bool readingInt = false;
        bool readingLen = false;
        bool readingStr = false;
        uint8_t byte;
        phr.seekg(start * sizeof(byte));
        while (phr.tellg() < end && phr.read((char *)(&byte), sizeof(byte)))
        {
            if (!readingInt && !readingLen && !readingStr)
            {
                switch (byte)
                {
                case 0x1a:
                    // String
                    readingStr = true;
                    readingLen = true;
                    break;
                case 0x02:
                    // int
                    readingInt = true;
                    readingLen = true;
                    break;
                }
            }
            else if (readingLen)
            {
                int len = 0;
                if (readingStr)
                {
                    if (byte & (1 << 7))
                    {
                        byte = (byte & ~(1 << 7)) | ((0 << 7) & (1 << 7));
                        phr.read((char *)&(len), sizeof(byte) * (int)byte);
                    }
                    else
                        len = (int)byte;
                    char buffer[len + 1];
                    phr.read((char *)&(buffer[0]), sizeof(byte) * len);
                    buffer[len] = 0;
                    readingStr = false;
                    if (title.empty())
                        title = buffer;
                    else
                        db = buffer;
                }
                else if (readingInt)
                {
                    len = (int)byte;
                    unsigned long intFound = 0;
                    uint8_t buffer = 0;
                    for (int j = 0; j < len; j++)
                    {
                        phr.read((char *)&buffer, sizeof(uint8_t));
                        intFound = intFound << 8 | (unsigned long)buffer;
                    }
                    readingInt = false;
                    if (id.empty())
                        id = to_string((int)intFound);
                    else
                        type = to_string((int)intFound);
                }
                readingLen = false;
            }
        }
        res = type + "|" + db + "|" + id + " " + title;
    }
    else
        cerr << "Impossible d'ouvrir le fichier" << endl;
}
