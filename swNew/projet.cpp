#include "Algorithm.h"
#include <chrono>

#ifdef _WIN32
#include <windows.h>
#elif MACOS
#include <sys/param.h>
#include <sys/sysctl.h>
#else
#include <unistd.h>
#endif

using namespace std;

int getNumCores();

int main(int argc, char *argv[])
{
	auto start = std::chrono::high_resolution_clock::now();
	if (argc < 3)
	{
		cout << "Error invalid argument" << endl
			 << "projet database query.fasta [-o IntOpenPenality] [-e IntExtensionPenality] [-b BLOSUM]" << endl
			 << "Default value : IntOpenPenality=11 IntExtensionPenality=1 BLOSUM=BLOSUM62" << endl;
		return EXIT_FAILURE;
	}

	clock_t before = clock();

	string pathBlosum = "BLOSUM62";
	int extension_gap = 1, open_gap = 11;
	int core = getNumCores();
	for (int i = 1; i < argc; ++i)
	{
		string current_arg = (string)argv[i];
		if (current_arg == "-o")
		{
			// Optional gap open penalty is set
			if (argv[i + 1] == NULL)
			{ 
				cout << "Invalid gap open penalty argument." << endl;
				exit(EXIT_FAILURE);
			}
			else
			{
				open_gap = atoi(argv[i + 1]);
				continue; // we can skip next argument as it is the gap open penalty value
			}
		}
		else if (current_arg == "-e")
		{
			// Optional gap expansion penalty is set
			current_arg = (string)argv[i + 1];
			if (argv[i + 1] == NULL)
			{ 
				cout << "Invalid gap expansion penalty argument." << endl;
				exit(EXIT_FAILURE);
			}
			else
			{
				extension_gap = atoi(argv[i + 1]);
				continue; // we can skip next argument as it is the gap expansion penalty value
			}
		}
		else if (current_arg == "-b")
		{
			// Optional BLOSUM matrix is set
			pathBlosum = (string)argv[i + 1];
			if (argv[i + 1] == NULL)
			{
				cout << "Invalid blosum argument." << endl;
				exit(EXIT_FAILURE);
			}
			continue;
		}
	}

	cout << "Number of CPU : " << core << endl;
	Algorithm *al = new Algorithm(argv[1], argv[2], pathBlosum, core, extension_gap, open_gap);
	//al->exactMatch();
	al->startMultithread();
	al->showResult();
	
	auto end = std::chrono::high_resolution_clock::now();
    auto sec_duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    auto ms_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nTime taken : " << sec_duration << "." << ms_duration % 1000  << "s" << std::endl;
	delete al;
	return 0;
}

int getNumCores()
{
#ifdef WIN32
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);
	return sysinfo.dwNumberOfProcessors;
#elif MACOS
	int nm[2];
	size_t len = 4;
	uint32_t count;

	nm[0] = CTL_HW;
	nm[1] = HW_AVAILCPU;
	sysctl(nm, 2, &count, &len, NULL, 0);

	if (count < 1)
	{
		nm[1] = HW_NCPU;
		sysctl(nm, 2, &count, &len, NULL, 0);
		if (count < 1)
		{
			count = 1;
		}
	}
	return count;
#else
	return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}

