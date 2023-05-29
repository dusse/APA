#include "ModelInitializer.hpp"



using namespace std;
using namespace chrono;


ModelInitializer::ModelInitializer(shared_ptr<Loader> load,
                                   shared_ptr<GridManager> grid):loader(move(load)), gridMng(move(grid)){
    logger.reset(new Logger());
    
    initialize();
    logger->writeMsg("[ModelInitializer] create...OK", DEBUG);
}

void ModelInitializer::initialize(){
    logger->writeMsg("[ModelInitializer] initialize...OK", DEBUG);
    
    int xRes = loader->resolution[0],
        yRes = loader->resolution[1],
        zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;

    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];

    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];

    double pfluence = 0.0;
    
    int i,j,k=1, idxG2;
    
    for( i=0; i<xRes+1; i++ ){
        for( j=0; j<yRes+1; j++ ){
                 
            idxG2 = IDX(i+1,j+1,k,xResG2,yResG2, zResG2);
                 
            pfluence = loader->getProtonFluence(i+1+int(domainShiftX/dx),j+1+int(domainShiftY/dy),k);
            double coord_x = domainShiftX + (i+1)*dx;
            double coord_y = domainShiftY + (j+1)*dy;

            double phi = (pow(coord_x,2) + pow(coord_y,2))/2;
            gridMng->setVectorVariableForNodeG2(idxG2, 
                VectorVar(PROTONFLUENCE, {coord_x, coord_y, pfluence, phi, 0.0, 0.0}));
            gridMng->setVectorVariableForNodeG2(idxG2, POTENTIAL, 0, phi);
            
        }
    }

    gridMng->sendBoundary2Neighbor(PROTONFLUENCE);
    gridMng->applyBC(PROTONFLUENCE);

    gridMng->sendBoundary2Neighbor(POTENTIAL);
    gridMng->applyBC(POTENTIAL);

        if (loader->neighbors2Send[NEIGHBOR_LEFT] == MPI_PROC_NULL){
            for( j = 0; j < yResG2; j++ ){
                double coord_y = domainShiftY + j*dy;
                double phi = (pow(coord_y,2))/2;
                idxG2 = IDX(0, j, k, xResG2, yResG2, zResG2);
                gridMng->setVectorVariableForNodeG2(idxG2, POTENTIAL, 0, phi);
                gridMng->setVectorVariableForNodeG2(idxG2, PROTONFLUENCE, 0, 0.0);
                gridMng->setVectorVariableForNodeG2(idxG2, PROTONFLUENCE, 1, coord_y);
            }
        }
        if (loader->neighbors2Send[NEIGHBOR_RIGHT] == MPI_PROC_NULL){
            for( j = 0; j < yResG2; j++ ){
                double coord_x = domainShiftX + (xRes+1)*dx;
                double coord_y = domainShiftY + j*dy;
                double phi = (pow(coord_x,2)+pow(coord_y,2))/2;
                idxG2 = IDX(xRes+1, j, k, xResG2, yResG2, zResG2);
                gridMng->setVectorVariableForNodeG2(idxG2, POTENTIAL, 0, phi);
                gridMng->setVectorVariableForNodeG2(idxG2, PROTONFLUENCE, 0, coord_x);
                gridMng->setVectorVariableForNodeG2(idxG2, PROTONFLUENCE, 1, coord_y);
            }
        }

        if (loader->neighbors2Send[NEIGHBOR_BOTTOM] == MPI_PROC_NULL){
            for( i = 0; i < xResG2; i++ ){
                double coord_x = domainShiftX + i*dx;
                double phi = (pow(coord_x,2))/2;
                idxG2 = IDX(i, 0, k, xResG2, yResG2, zResG2);
                gridMng->setVectorVariableForNodeG2(idxG2, POTENTIAL, 0, phi);
                gridMng->setVectorVariableForNodeG2(idxG2, PROTONFLUENCE, 0, coord_x);
                gridMng->setVectorVariableForNodeG2(idxG2, PROTONFLUENCE, 1, 0.0);
            }
        }

        if (loader->neighbors2Send[NEIGHBOR_TOP] == MPI_PROC_NULL){
            for( i = 0; i < xResG2; i++ ){
                double coord_x = domainShiftX + i*dx;
                double coord_y = domainShiftY + (yRes+1)*dy;
                double phi = (pow(coord_x,2)+pow(coord_y,2))/2;
                idxG2 = IDX(i, yRes+1, k, xResG2, yResG2, zResG2);
                gridMng->setVectorVariableForNodeG2(idxG2, POTENTIAL, 0, phi);
                gridMng->setVectorVariableForNodeG2(idxG2, PROTONFLUENCE, 0, coord_x);
                gridMng->setVectorVariableForNodeG2(idxG2, PROTONFLUENCE, 1, coord_y);
            }
        }
    // gridMng->applyBC(POTENTIAL);
}

                          
ModelInitializer::~ModelInitializer(){
    finilize();
    logger->writeMsg("[ModelInitializer] delete...OK", DEBUG);
}

void ModelInitializer::finilize(){
    
}
