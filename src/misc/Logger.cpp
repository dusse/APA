#include "Logger.hpp"

#include <iostream>
#include <string>
using namespace std;



Logger::Logger()
{
}

void Logger::writeMsg(const char* input, int level)
{

    int rank ;
    
    set<int> allowedRanks = {0};
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (level <= MINIMAL_LEVEL && allowedRanks.count(rank)){
        cout << "[" << rank << "] " << "[" << level << "] " << input <<  "." << endl;
    }
    
    if (level == CRITICAL){
        cout << "[" << rank << "] " << "[" << level << "] " << input <<  "." << endl;
    }
}

void Logger::writeMsg(const char* input)
{
    int rank ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0){
        cout << "[INFO] " << input <<  "." << endl;
    }
}



Logger::~Logger()
{
}


