#include "GridManager.hpp"

using namespace std;
using namespace chrono;

GridManager::GridManager(shared_ptr<Loader> loader){
    this->loader=loader;
    logger.reset(new Logger());
    initialize();
    string msg ="[GridManager] init...OK";
    logger->writeMsg(msg.c_str(), DEBUG);
}


GridManager::~GridManager(){
    delete nodesG2vars;
    
    for (int t = 0; t < 27; t++) {
        delete [] sendIdx[t];
        delete [] recvIdx[t];
    }
}


void GridManager::initialize(){
    logger->writeMsg("[GridManager] initialize() ...", DEBUG);
    
    int xRes = loader->resolution[0];
    int yRes = loader->resolution[1];
    int zRes = loader->resolution[2];
    
    this->G2nodesNumber = (xRes+2)*(yRes+2)*(zRes+2);
    
    nodesG2vars = new VectorVar*[G2nodesNumber*SIZEG2];

    initG2Nodes();    

    sendRecvIndecis4MPI();
    sendRecvIndecis4MPIonG4();
    
    initBoundaryIndecies();

    getVectorVariablesForAllNodes();
    logger->writeMsg("[GridManager] initialize() ...OK", DEBUG);
}


void GridManager::initG2Nodes(){
    auto start_time = high_resolution_clock::now();
    logger->writeMsg("[GridManager] init Node start ", DEBUG);
    int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    
    int i,j,k,idx;
    for( i=0; i<xResG2; i++ ){
        for( j=0; j<yResG2; j++ ){
            for( k=0; k<zResG2; k++ ){
                idx = IDX(i,j,k,xResG2,yResG2,zResG2);
                VectorVar* fluence = new VectorVar(PROTONFLUENCE, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
                nodesG2vars[G2nodesNumber*PROTONFLUENCE+idx] = fluence;
                VectorVar* pot = new VectorVar(POTENTIAL, {0.0});
                nodesG2vars[G2nodesNumber*POTENTIAL+idx] = pot;                
                VectorVar* outvar = new VectorVar(OUTPUTVAR, {0.0, 0.0});
                nodesG2vars[G2nodesNumber*OUTPUTVAR+idx] = outvar;
                VectorVar* auxvar1 = new VectorVar(AUXVAR1, {0.0, 0.0});
                nodesG2vars[G2nodesNumber*AUXVAR1+idx] = auxvar1;
                VectorVar* auxvar2 = new VectorVar(AUXVAR2, {0.0, 0.0, 0.0});
                nodesG2vars[G2nodesNumber*AUXVAR2+idx] = auxvar2;
            }
        }
    }
    auto end_time = high_resolution_clock::now();
    string msg ="[GridManager] init Node G2 duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}

void GridManager::initBoundaryIndecies(){
    int xRes = loader->resolution[0];
    int yRes = loader->resolution[1];
    int zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    int i,j,k, idx2set, idx2use;

    if( xRes != 1 ){
        
        if (loader->neighbors2Send[NEIGHBOR_LEFT] == MPI_PROC_NULL){
            for ( j=0; j<yResG2; j++){
                for ( k=0; k<zResG2; k++){
                    idx2use = IDX(1,j,k, xResG2,yResG2,zResG2);
                    idx2set = IDX(0,j,k, xResG2,yResG2,zResG2);
                        
                    idxs4BoundaryX2fill[idx2use] = idx2set;
                }
            }
        }
        
        if (loader->neighbors2Send[NEIGHBOR_RIGHT] == MPI_PROC_NULL){
            for ( j=0; j<yResG2; j++){
                for ( k=0; k<zResG2; k++){
                    idx2use = IDX(xResG2-2,j,k, xResG2,yResG2,zResG2);
                    idx2set = IDX(xResG2-1,j,k, xResG2,yResG2,zResG2);
                    
                    idxs4BoundaryX2fill[idx2use] = idx2set;
                }
            }
        }
    }
    
    if( yRes != 1 ){
        
        if (loader->neighbors2Send[NEIGHBOR_BOTTOM] == MPI_PROC_NULL){
            for ( i=0; i<xResG2; i++){
                for ( k=0; k<zResG2; k++){
                    idx2use = IDX(i,1,k, xResG2,yResG2,zResG2);
                    idx2set = IDX(i,0,k, xResG2,yResG2,zResG2);
                    
                    idxs4BoundaryY2fill[idx2use] = idx2set;
                }
            }
        }
        if (loader->neighbors2Send[NEIGHBOR_TOP] == MPI_PROC_NULL){
            for ( i=0; i<xResG2; i++){
                for ( k=0; k<zResG2; k++){
                    idx2use = IDX(i,yResG2-2,k, xResG2,yResG2,zResG2);
                    idx2set = IDX(i,yResG2-1,k, xResG2,yResG2,zResG2);
                    
                    idxs4BoundaryY2fill[idx2use] = idx2set;
                }
            }
        }
    }
    
    if( zRes != 1 ){
        
        if (loader->neighbors2Send[NEIGHBOR_BACK] == MPI_PROC_NULL){
            for ( i=0; i<xResG2; i++){
                for ( j=0; j<yResG2; j++){
                    idx2use = IDX(i,j,1, xResG2,yResG2,zResG2);
                    idx2set = IDX(i,j,0, xResG2,yResG2,zResG2);
                    
                    idxs4BoundaryZ2fill[idx2use] = idx2set;
                }
            }
        }
        if (loader->neighbors2Send[NEIGHBOR_FRONT] == MPI_PROC_NULL){
            for ( i=0; i<xResG2; i++){
                for ( j=0; j<yResG2; j++){
                    idx2use = IDX(i,j,zResG2-2, xResG2,yResG2,zResG2);
                    idx2set = IDX(i,j,zResG2-1, xResG2,yResG2,zResG2);
                    
                    idxs4BoundaryZ2fill[idx2use] = idx2set;
                }
            }
        }
    }
    
    string msg ="[GridManager] initialized idxs4BoundaryX2fill.size() = "
    +to_string(idxs4BoundaryX2fill.size())
    +"\n                       idxs4BoundaryY2fill.size() = "
    +to_string(idxs4BoundaryY2fill.size())
    +"\n                       idxs4BoundaryZ2fill.size() = "
    +to_string(idxs4BoundaryZ2fill.size());
    logger->writeMsg(msg.c_str(), DEBUG);
}

void GridManager::applyBC(int varName){
    //nothing for periodic BC, maps are empty
    auto start_time = high_resolution_clock::now();
    
    int varDim = nodesG2vars[G2nodesNumber*varName]->getSize();
    int varShift = G2nodesNumber*varName;
    int dim;
    const double * vectorVar;
    
    for ( const auto &keyval : idxs4BoundaryX2fill ) {
        vectorVar = nodesG2vars[varShift+keyval.first]->getValue();
        for ( dim = 0; dim < varDim; dim++ ){
            nodesG2vars[varShift+keyval.second]->setValue( dim, vectorVar[dim]);
        }
    }
    
    for ( const auto &keyval : idxs4BoundaryY2fill ) {
        vectorVar = nodesG2vars[varShift+keyval.first]->getValue();
        for ( dim = 0; dim < varDim; dim++ ){
            nodesG2vars[varShift+keyval.second]->setValue( dim, vectorVar[dim]);
        }
    }
    
    for ( const auto &keyval : idxs4BoundaryZ2fill ) {
        vectorVar = nodesG2vars[varShift+keyval.first]->getValue();
        for ( dim = 0; dim < varDim; dim++ ){
            nodesG2vars[varShift+keyval.second]->setValue( dim, vectorVar[dim]);
        }
    }
    
    auto end_time = high_resolution_clock::now();
    string msg ="[GridManager] apply the same value BC for var = "+to_string(varName) +"  duration = "
    +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
    
}

void GridManager::sendRecvIndecis4MPI(){
    int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    int idx, i, j, k, a, b, c, t;
    vector<int> idxX, idxY, idxZ;
    
    for (a = -1; a <= 1; a++) {
        for (b = -1; b <= 1; b++) {
            for (c = -1; c <= 1; c++) {
                t = (1 + c) + 3 * ((1 + b) + 3 * (1 + a));

                switch (a) {
                    case -1: idxX = {1   , 1   , xRes + 1, xRes + 1}; break;
                    case  0: idxX = {1   , xRes, 1       , xRes    }; break;
                    case  1: idxX = {xRes, xRes, 0       , 0       }; break;
                }
                switch (b) {
                    case -1: idxY = {1   , 1   , yRes + 1, yRes + 1}; break;
                    case  0: idxY = {1   , yRes, 1       , yRes    }; break;
                    case  1: idxY = {yRes, yRes, 0       , 0       }; break;
                }
                switch (c) {
                    case -1: idxZ = {1   , 1   , zRes + 1, zRes + 1}; break;
                    case  0: idxZ = {1   , zRes, 1       , zRes    }; break;
                    case  1: idxZ = {zRes, zRes, 0       , 0       }; break;
                }
                
                counter[t] = (idxX[1] - idxX[0] + 1)*
                             (idxY[1] - idxY[0] + 1)*
                             (idxZ[1] - idxZ[0] + 1);
                sendIdx[t] = new int[counter[t]];
                recvIdx[t] = new int[counter[t]];
                
                int lineIDX = 0;
                for( i = idxX[0]; i <= idxX[1]; i++ ) {
                    for( j = idxY[0]; j <= idxY[1]; j++ ) {
                        for( k = idxZ[0]; k <= idxZ[1]; k++ ) {
                            idx = IDX(i,j,k, xResG2, yResG2, zResG2);
                            sendIdx[t][lineIDX] = idx;
                            lineIDX++;
                        }
                    }
                }

                lineIDX = 0;
                for (i = idxX[2]; i <= idxX[3]; i++) {
                    for (j = idxY[2]; j <= idxY[3]; j++) {
                        for (k = idxZ[2]; k <= idxZ[3]; k++) {
                            idx = IDX(i,j,k, xResG2, yResG2, zResG2);
                            recvIdx[t][lineIDX] = idx;
                            lineIDX++;
                        }
                    }
                }
            }
        }
    }
}


void GridManager::sendBoundary2Neighbor(int varName){
    
    auto start_time = high_resolution_clock::now();
    int varShift = G2nodesNumber*varName;
    int varDim = nodesG2vars[varShift]->getSize();
    int t, i;
    double *sendBuf[27];
    double *recvBuf[27];
    const double* vectorVar;
    
    for( t = 0; t < 27; t++ ){
        sendBuf[t] = new double[counter[t]*varDim*sizeof(double)];
        recvBuf[t] = new double[counter[t]*varDim*sizeof(double)];
        for( i = 0; i < counter[t]; i++ ){
            vectorVar = nodesG2vars[varShift+sendIdx[t][i]]->getValue();
            for( int dim = 0; dim < varDim; dim++ ){
                sendBuf[t][varDim*i+dim] = vectorVar[dim];
            }
        }
    }
    
    MPI_Status st;
    for( t = 0; t < 27; t++ ){
        if( t != 13 ){
            MPI_Sendrecv(sendBuf[t], counter[t]*varDim, MPI_DOUBLE, loader->neighbors2Send[t], t,
                         recvBuf[t], counter[t]*varDim, MPI_DOUBLE, loader->neighbors2Recv[t], t,
                         MPI_COMM_WORLD, &st);
        }
    }
    
    for( t = 0; t < 27; t++ ){
        if( t != 13 && loader->neighbors2Recv[t] != MPI_PROC_NULL ){
            for( i = 0; i < counter[t]; i++ ){
                for( int dim = 0; dim < varDim; dim++ ){
                    nodesG2vars[varShift+recvIdx[t][i]]
                    ->setValue( dim, recvBuf[t][varDim*i+dim]);
                }
            }
        }
    }
    
    
    for( t = 0; t < 27; t++ ){
        delete [] sendBuf[t];
        delete [] recvBuf[t];
    }
    
    auto end_time = high_resolution_clock::now();
    string msg ="[GridManager] send boundary for var = "+to_string(varName)
                +" duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}




void GridManager::sendRecvIndecis4MPIonG4(){
    int xRes = loader->resolution[0],
    yRes = loader->resolution[1],
    zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    int xResG4 = xRes+4, yResG4 = yRes+4, zResG4 = zRes+4;
    int idx, i, j, k, a, b, c, t;
    vector<int> idxX, idxY, idxZ;
    
    for (a = -1; a <= 1; a++) {
        for (b = -1; b <= 1; b++) {
            for (c = -1; c <= 1; c++) {
                
                switch (a) {
                    case -1: idxX = {3   , 3       , xRes + 3, xRes + 3}; break;
                    case  0: idxX = {1   , xRes + 2, 1       , xRes + 2}; break;
                    case  1: idxX = {xRes, xRes    , 0       , 0       }; break;
                }
                switch (b) {
                    case -1: idxY = {3   , 3       , yRes + 3, yRes + 3}; break;
                    case  0: idxY = {1   , yRes + 2, 1       , yRes + 2}; break;
                    case  1: idxY = {yRes, yRes    , 0       , 0       }; break;
                }
                switch (c) {
                    case -1: idxZ = {3   , 3       , zRes + 3, zRes + 3}; break;
                    case  0: idxZ = {1   , zRes + 2, 1       , zRes + 2}; break;
                    case  1: idxZ = {zRes, zRes    , 0       , 0       }; break;
                }
                
                t = (1 + c) + 3 * ((1 + b) + 3 * (1 + a));
                
                counterOnG4[t] = (idxX[1] - idxX[0] + 1)*
                (idxY[1] - idxY[0] + 1)*
                (idxZ[1] - idxZ[0] + 1);
                sendIdxOnG4[t] = new int[counterOnG4[t]];
                recvIdxOnG4[t] = new int[counterOnG4[t]];
                
                int lineIDX = 0;
                for (i = idxX[0]; i <= idxX[1]; i++) {
                    for (j = idxY[0]; j <= idxY[1]; j++) {
                        for (k = idxZ[0]; k <= idxZ[1]; k++) {
                            idx = IDX(i,j,k, xResG4, yResG4, zResG4);
                            sendIdxOnG4[t][lineIDX] = idx;
                            lineIDX++;
                        }
                    }
                }
                lineIDX = 0;
                for (i = idxX[2]; i <= idxX[3]; i++) {
                    for (j = idxY[2]; j <= idxY[3]; j++) {
                        for (k = idxZ[2]; k <= idxZ[3]; k++) {
                            idx = IDX(i,j,k, xResG4, yResG4, zResG4);
                            recvIdxOnG4[t][lineIDX] = idx;
                            lineIDX++;
                        }
                    }
                }
            }
        }
    }
}



void GridManager::smooth(int varName){
    
    int i, j, k, idx, idxG4;
    int xRes = loader->resolution[0],
    yRes = loader->resolution[1],
    zRes = loader->resolution[2];
    int xResG4 = xRes+4, yResG4 = yRes+4, zResG4 = zRes+4;
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    int totG4 = xResG4*yResG4*zResG4;
    auto start_time = high_resolution_clock::now();
    int varDim = nodesG2vars[G2nodesNumber*varName]->getSize();
    int t;
    double *sendBuf[27];
    double *recvBuf[27];
    const double* vectorVar;
    
    double* varValues = new double[totG4*varDim*sizeof(double)];
    
    for ( i=0; i<xResG2; i++){
        for ( j=0; j<yResG2; j++){
            for ( k=0; k<zResG2; k++){
                idx   = IDX(i  ,j  ,k  ,xResG2,yResG2,zResG2);
                idxG4 = IDX(i+1,j+1,k+1,xResG4,yResG4,zResG4);
                vectorVar = nodesG2vars[G2nodesNumber*varName+idx]->getValue();
                for (int dim=0;dim<varDim;dim++) {
                    varValues[varDim*idxG4+dim] = vectorVar[dim];
                }
            }
        }
    }
    
    for (t = 0; t < 27; t++) {
        sendBuf[t] = new double[counterOnG4[t]*varDim*sizeof(double)];
        recvBuf[t] = new double[counterOnG4[t]*varDim*sizeof(double)];
        for(i = 0; i < counterOnG4[t]; i++){
            int si = sendIdxOnG4[t][i];
            for (int dim=0;dim<varDim;dim++) {
                sendBuf[t][varDim*i+dim] = varValues[varDim*si+dim];
                
            }
        }
    }
    
    MPI_Status st;
    for (t = 0; t < 27; t++) {
        if (t != 13){
            MPI_Sendrecv(sendBuf[t], counterOnG4[t]*varDim, MPI_DOUBLE, loader->neighbors2Send[t], t,
                         recvBuf[t], counterOnG4[t]*varDim, MPI_DOUBLE, loader->neighbors2Recv[t], t,
                         MPI_COMM_WORLD, &st);
        }
    }
    for (t = 0; t < 27; t++) {
        if (t != 13 && loader->neighbors2Recv[t] != MPI_PROC_NULL) {
            for(i = 0; i < counterOnG4[t]; i++){
                int ri = recvIdxOnG4[t][i];
                for (int dim=0;dim<varDim;dim++) {
                    varValues[varDim*ri+dim] = recvBuf[t][varDim*i+dim];
                }
            }
        }
    }
    
    
    for (t = 0; t < 27; t++) {
        delete [] sendBuf[t];
        delete [] recvBuf[t];
    }
    
    int zeroOrderNeighb[6][3]  =
    {{-1,0 ,0 }, {+1,0 ,0 },
        {0 ,-1,0 }, {0 ,+1,0 },
        {0 ,0 ,-1}, {0 ,0 ,+1}};
    
    int firstOrderNeighb[12][3] =
    {{-1,-1, 0}, {-1,+1, 0},
        {-1, 0,-1}, {-1, 0,+1},
        {0 ,-1,-1}, { 0,-1,+1},
        {0 ,+1,-1}, { 0,+1,+1},
        {+1, 0,-1}, {+1, 0,+1},
        {+1,-1, 0}, {+1,+1, 0}};
    
    int secndOrderNeighb[8][3] =
    {{-1,-1,-1}, {-1,-1,+1},
        {-1,+1,-1}, {-1,+1,+1},
        {+1,-1,-1}, {+1,-1,+1},
        {+1,+1,-1}, {+1,+1,+1}};
    
    const double k2 = 0.125;// 1/8
    const double k3 = 0.0625;// 1/16
    const double k4 = 0.03125;// 1/32
    const double k5 = 0.015625;// 1/64
    //                    +--------
    //  0.125*1+0.0625*6+0.03125*12+0.015625*8 = 1
    
    int foi, soi, toi;
    int idxFo, idxSo, idxTo;
    
    for ( i=1; i<xResG2+1; i++){
        for ( j=1; j<yResG2+1; j++){
            for ( k=1; k<zResG2+1; k++){
                idxG4 = IDX(i,j,k,xResG4,yResG4,zResG4);
                
                for (int dim=0;dim<varDim;dim++) {
                    
                    double smoothedVal = k2*varValues[varDim*idxG4+dim];
                    
                    for(foi = 0; foi<6; foi++){
                        idxFo = IDX(i+zeroOrderNeighb[foi][0],
                                    j+zeroOrderNeighb[foi][1],
                                    k+zeroOrderNeighb[foi][2],
                                    xResG4,yResG4,zResG4);
                        smoothedVal += k3*varValues[varDim*idxFo+dim];
                    }
                    
                    for(soi = 0; soi<12; soi++){
                        idxSo = IDX(i+firstOrderNeighb[soi][0],
                                    j+firstOrderNeighb[soi][1],
                                    k+firstOrderNeighb[soi][2],
                                    xResG4,yResG4,zResG4);
                        smoothedVal += k4*varValues[varDim*idxSo+dim];
                    }
                    
                    for(toi = 0; toi<8; toi++){
                        idxTo = IDX(i+secndOrderNeighb[toi][0],
                                    j+secndOrderNeighb[toi][1],
                                    k+secndOrderNeighb[toi][2],
                                    xResG4,yResG4,zResG4);
                        smoothedVal += k5*varValues[varDim*idxTo+dim];
                    }
                    idx = IDX(i-1,j-1,k-1,xResG2,yResG2,zResG2);
                    nodesG2vars[G2nodesNumber*varName+idx]->setValue(dim, smoothedVal);
                }
            }
        }
    }
    
    
    delete varValues;
    auto end_time = high_resolution_clock::now();
    string msg ="[GridManager] smooth "+to_string(varName)
    +" duration = "+to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
    
}



vector<vector<VectorVar>> GridManager::getVectorVariablesForAllNodes(){
    
    int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    vector<vector<VectorVar>> result;
    result.reserve(xRes*yRes*zRes);
    int i,j,k,idxG2;
    
    
    for ( i=1; i<xRes+1; i++){
        for ( j=1; j<yRes+1; j++){
            for ( k=1; k<zRes+1; k++){
                idxG2 = IDX(i,j,k,xResG2,yResG2,zResG2);
                vector<VectorVar> allVars;
                allVars.push_back(*nodesG2vars[G2nodesNumber*PROTONFLUENCE+idxG2]);                
                allVars.push_back(*nodesG2vars[G2nodesNumber*OUTPUTVAR+idxG2]);                
                result.push_back(allVars);
            }
        }
    }
    
    return result;
}



void GridManager::setVectorVariableForNodeG2(int idx, VectorVar variable){
    nodesG2vars[G2nodesNumber*variable.getName()+idx]->setValue(variable.getValue());
}

void GridManager::setVectorVariableForNodeG2(int idx, int name, int dim, double value){
    nodesG2vars[G2nodesNumber*name+idx]->setValue(dim, value);
}


void GridManager::addVectorVariableForNodeG2(int idx, int name, int dim, double value){
    nodesG2vars[G2nodesNumber*name+idx]->addValue(dim, value);
}


VectorVar** GridManager::getVectorVariableOnG2(int varName){
     return &nodesG2vars[G2nodesNumber*varName];
}

