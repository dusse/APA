#ifndef GridManager_hpp
#define GridManager_hpp
#include <stdio.h>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <memory>
#include <chrono>
#include <mpi.h>
#include "../misc/Logger.hpp"
#include "../misc/Misc.hpp"
#include "../input/Loader.hpp"

#include "../common/variables/VectorVar.hpp"


enum G2VAR{
    PROTONFLUENCE,
    POTENTIAL,
    OUTPUTVAR,
    AUXVAR1,
    AUXVAR2,
    SIZEG2
};

//# for boundary conditions
#define NEIGHBOR_LEFT   4
#define NEIGHBOR_RIGHT  22
#define NEIGHBOR_BOTTOM 10
#define NEIGHBOR_TOP    16
#define NEIGHBOR_BACK   12
#define NEIGHBOR_FRONT  14

class GridManager{
    
    
private:
    
    std::unique_ptr<Logger> logger;
    std::shared_ptr<Loader> loader;
    
    std::map<int, int> idxs4BoundaryX2fill;
    std::map<int, int> idxs4BoundaryY2fill;
    std::map<int, int> idxs4BoundaryZ2fill;
    
    int G2nodesNumber;
    
    VectorVar** nodesG2vars;
    
    int counter[27], counterOnG4[27];
    int* sendIdx[27];
    int* recvIdx[27];
    int* sendIdxOnG4[27];
    int* recvIdxOnG4[27];
    
    void initialize();
    void initG2Nodes();
    void initBoundaryIndecies();
    
    void sendRecvIndecis4MPI();
    void sendRecvIndecis4MPIonG4();
    
public:
    
    GridManager(std::shared_ptr<Loader>);
    ~GridManager();
    
    VectorVar** getVectorVariableOnG2(int);
        
    std::vector<std::vector<VectorVar>> getVectorVariablesForAllNodes();
    
    void setVectorVariableForNodeG2(int, VectorVar);
    
    void setVectorVariableForNode(int, VectorVar);
    void setVectorVariableForNodeG2(int , int , int, double);
    void addVectorVariableForNodeG2(int , int , int, double);
    
    void sendBoundary2Neighbor(int);
    void smooth(int);
    void applyBC(int);
    
    
};
#endif /* GridManager_hpp */
