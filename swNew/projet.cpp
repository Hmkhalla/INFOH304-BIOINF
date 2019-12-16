#include "Algorithm.h"
//#include <pthread.h>

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
	int core = getNumCores();
	cout << "Number of cores : "<<core<<endl;
	Algorithm *al = new Algorithm(argv[1],argv[2],argv[3], core);
	//al->exactMatch();
	//al->swAlgo();
	cout << "1"<< endl;
	al->startMultithread();

	delete al;
	return 0;
}

int getNumCores() {
#ifdef WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    return sysinfo.dwNumberOfProcessors;
#elif MACOS
    int nm[2];
    size_t len = 4;
    uint32_t count;

    nm[0] = CTL_HW; nm[1] = HW_AVAILCPU;
    sysctl(nm, 2, &count, &len, NULL, 0);

    if(count < 1) {
        nm[1] = HW_NCPU;
        sysctl(nm, 2, &count, &len, NULL, 0);
        if(count < 1) { count = 1; }
    }
    return count;
#else
    return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}
