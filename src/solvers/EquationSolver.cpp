#include "EquationSolver.hpp"



using namespace std;
using namespace chrono;




EquationSolver::EquationSolver(shared_ptr<Loader> load,
                           shared_ptr<GridManager> grid):loader(move(load)),
                            gridMng(move(grid)){
    logger.reset(new Logger());
    
    initDampingCoeff();
    initialize();    

    logger->writeMsg("[EquationSolver] create...OK", DEBUG);
}

void EquationSolver::initDampingCoeff(){

    int xRes = loader->resolution[0],
        yRes = loader->resolution[1],
        zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    
    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];
    
    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];
    
    double Lx = loader->boxSizes[0]/dx,
           Ly = loader->boxSizes[1]/dy,
           Lz = loader->boxSizes[2]/dz;
    
    int idx2use;
    int nG2 = xResG2*yResG2*zResG2;
    int i, j, k=1, idx;
    
    dampingCoeff = new double[nG2*sizeof(double)];
    
    for( idx = 0; idx < nG2; idx++ ){
        dampingCoeff[idx] = 1.0;
    }
    
    double leftXWidth = loader->dampingBoundaryWidth[0][0];
    double rigtXWidth = loader->dampingBoundaryWidth[0][1];
    
    double leftYWidth = loader->dampingBoundaryWidth[1][0];
    double rigtYWidth = loader->dampingBoundaryWidth[1][1];
    
    double normLengthXleft, normLengthXrigt,
           normLengthYleft, normLengthYrigt, minCoef;
    
    double x, y;
    
    for( i = 0; i < xResG2; i++ ){
        for( j = 0; j < yResG2; j++ ){
                
                x = i + domainShiftX/dx;
                y = j + domainShiftY/dy;
                
                idx2use = IDX(i,j,k, xResG2,yResG2,zResG2);
                normLengthXleft = 1; normLengthXrigt = 1;
                normLengthYleft = 1; normLengthYrigt = 1;
                minCoef = 1;
                
                
                if( leftXWidth > 0.0 ){
                    normLengthXleft = x/leftXWidth*exp(x-leftXWidth);
                    if( normLengthXleft >= 0 && normLengthXleft < minCoef ){
                        minCoef = normLengthXleft;
                    }
                }
                
                if( rigtXWidth > 0.0 ){
                    normLengthXrigt = ( Lx + dx - x )/rigtXWidth*exp(Lx + dx - x - rigtXWidth);
                    if( normLengthXrigt >= 0 && normLengthXrigt < minCoef ){
                        minCoef = normLengthXrigt;
                    }
                }
                
                if( leftYWidth > 0.0 ){
                    normLengthYleft = y/leftYWidth*exp(y-leftYWidth);
                    if( normLengthYleft >= 0 && normLengthYleft < minCoef){
                        minCoef = normLengthYleft;
                    }
                }
                
                if( rigtYWidth > 0.0 ){
                    normLengthYrigt = ( Ly + dy - y )/rigtYWidth*exp(Ly + dy - y - rigtYWidth);
                    if( normLengthYrigt >= 0 && normLengthYrigt < minCoef ){
                        minCoef = normLengthYrigt;
                    }
                }
                
                dampingCoeff[idx2use] = minCoef;
            
        }
    }
}

void EquationSolver::initialize(){
    int xTotSize = loader->totPixelsPerBoxSide[0]+1;
    int yTotSize = loader->totPixelsPerBoxSide[1]+1;
    
    initialProtonFluence = new double[xTotSize*yTotSize*sizeof(double)];
    bkgrndProtonFluence  = new double[xTotSize*yTotSize*sizeof(double)];

    int i,j,k=1, idxGlob;
    
    for( i=0; i<xTotSize; i++ ){
        for( j=0; j<yTotSize; j++ ){                 
            idxGlob = j + xTotSize*i;           
            double pfluence = loader->getProtonFluence(i,j,k);
            initialProtonFluence[idxGlob] = pfluence;
            pfluence = loader->getBackgroundProtonFluence(i,j,k); 
            if( pfluence <= 0.0 ){
                throw runtime_error("check input file!!! imposible background proton fluence value = "+to_string(pfluence));                
            }
            bkgrndProtonFluence[idxGlob] = pfluence;
        }
    }

    solve(-1);
    logger->writeMsg("[EquationSolver] initialize...OK", DEBUG);
}


int EquationSolver::solve(double time){

    int rank ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int xRes = loader->resolution[0],
        yRes = loader->resolution[1],
        zRes = loader->resolution[2];
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    int nG2= xResG2*yResG2*zResG2;

    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];

    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];

    string msg ="[EquationSolver] start dx ="+to_string(dx)+"; dy = "+to_string(dy);
    logger->writeMsg(msg.c_str(), DEBUG);

    msg ="[EquationSolver]  domainShiftX ="+to_string(domainShiftX)
    +"; domainShiftY = "+to_string(domainShiftY);
    logger->writeMsg(msg.c_str(), DEBUG);

    VectorVar** pfluence = gridMng->getVectorVariableOnG2(PROTONFLUENCE);
    VectorVar** pot = gridMng->getVectorVariableOnG2(POTENTIAL);
    VectorVar** aux1 = gridMng->getVectorVariableOnG2(AUXVAR1);
    VectorVar** aux2 = gridMng->getVectorVariableOnG2(AUXVAR2);
    double* phi = new double[nG2*sizeof(double)];
    double* phi_dx = new double[nG2*sizeof(double)];
    double* phi_dy = new double[nG2*sizeof(double)];
    double* phi_dxdx = new double[nG2*sizeof(double)];
    double* phi_dxdy = new double[nG2*sizeof(double)];
    double* phi_dydy = new double[nG2*sizeof(double)];
    double* det = new double[nG2*sizeof(double)];

    for( int ijkG2 = 0; ijkG2 < nG2; ijkG2++ ){
        phi[ijkG2] = pot[ijkG2]->getValue()[0];
        phi_dx[ijkG2] = 0.0;
        phi_dy[ijkG2] = 0.0;
        phi_dxdx[ijkG2] = 0.0;
        phi_dxdy[ijkG2] = 0.0;
        phi_dydy[ijkG2] = 0.0;
    }
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int idxXleft, idxXrigt, idxYleft, idxYrigt;
    int i,j,k=1,idxOnG2;
    try{

        
        for( i = 1; i < xRes + 1; i++ ){
            for( j = 1; j < yRes + 1; j++ ){

                    idxOnG2  = IDX(i, j, k, xResG2, yResG2, zResG2);

                    idxXleft = IDX(i-1, j, k, xResG2, yResG2, zResG2);
                    idxXrigt = IDX(i+1, j, k, xResG2, yResG2, zResG2);
                    phi_dx[idxOnG2] = 0.5*(phi[idxXrigt]-phi[idxXleft])/dx;

                    idxYleft = IDX(i, j-1, k, xResG2, yResG2, zResG2);
                    idxYrigt = IDX(i, j+1, k, xResG2, yResG2, zResG2);
                    phi_dy[idxOnG2] = 0.5*(phi[idxYrigt]-phi[idxYleft])/dy;
            }
        }
        int idxRight, idxLeft;

        for( int ijkG2 = 0; ijkG2 < nG2; ijkG2++ ){
            gridMng->setVectorVariableForNodeG2(ijkG2, AUXVAR1, 0, phi_dx[ijkG2]);
            gridMng->setVectorVariableForNodeG2(ijkG2, AUXVAR1, 1, phi_dy[ijkG2]);                       
         }

        gridMng->sendBoundary2Neighbor(AUXVAR1);

        for( int ijkG2 = 0; ijkG2 < nG2; ijkG2++ ){
            phi_dx[ijkG2] = aux1[ijkG2]->getValue()[0];
            phi_dy[ijkG2] = aux1[ijkG2]->getValue()[1];
        }


        //****** serve box boundaries ******////
        if (loader->neighbors2Send[NEIGHBOR_LEFT] == MPI_PROC_NULL){

            for( j = 0; j < yResG2; j++ ){
                idxOnG2  = IDX(0, j, k, xResG2, yResG2, zResG2);
                //this is so because of the initial potential form
                phi_dx[idxOnG2] = pfluence[idxOnG2]->getValue()[0];

                if ( j > 0 && j < (yResG2-1) ){
                    idxYleft = IDX(0, j-1, k, xResG2, yResG2, zResG2);
                    idxYrigt = IDX(0, j+1, k, xResG2, yResG2, zResG2);
                    phi_dy[idxOnG2] = 0.5*(phi[idxYrigt]-phi[idxYleft])/dy;
                }else{
                    phi_dy[idxOnG2] = pfluence[idxOnG2]->getValue()[1];
                }
            }
        }


        if (loader->neighbors2Send[NEIGHBOR_RIGHT] == MPI_PROC_NULL){
            for( j = 0; j < yResG2; j++ ){
                idxOnG2 = IDX(xRes+1, j, k, xResG2, yResG2, zResG2);
                int idxOnG2close  = IDX(xRes, j, k, xResG2, yResG2, zResG2);
                phi_dx[idxOnG2] = pfluence[idxOnG2]->getValue()[0];

                if ( j > 0 && j < (yResG2-1) ){
                    idxYleft = IDX(xRes+1, j-1, k, xResG2, yResG2, zResG2);
                    idxYrigt = IDX(xRes+1, j+1, k, xResG2, yResG2, zResG2);
                    phi_dy[idxOnG2] = 0.5*(phi[idxYrigt]-phi[idxYleft])/dy;
                }else{
                    phi_dy[idxOnG2] = pfluence[idxOnG2]->getValue()[1];
                }
            }
        }

        if (loader->neighbors2Send[NEIGHBOR_BOTTOM] == MPI_PROC_NULL){
            for( i = 0; i < xResG2; i++ ){
                idxOnG2  = IDX(i, 0, k, xResG2, yResG2, zResG2);
                if ( i > 0 && i < (xResG2-1) ){
                    idxXleft = IDX(i-1, 0, k, xResG2, yResG2, zResG2);
                    idxXrigt = IDX(i+1, 0, k, xResG2, yResG2, zResG2);
                    phi_dx[idxOnG2] = 0.5*(phi[idxXrigt]-phi[idxXleft])/dx;
                }else{
                    phi_dx[idxOnG2] = pfluence[idxOnG2]->getValue()[0];
                }
                phi_dy[idxOnG2] = pfluence[idxOnG2]->getValue()[1];
            }
        }

        if (loader->neighbors2Send[NEIGHBOR_TOP] == MPI_PROC_NULL){
            for( i = 0; i < xResG2; i++ ){
                idxOnG2 = IDX(i, yRes+1, k, xResG2, yResG2, zResG2);
                if ( i > 0 && i < (xResG2-1) ){
                    idxXleft = IDX(i-1, yRes+1, k, xResG2, yResG2, zResG2);
                    idxXrigt = IDX(i+1, yRes+1, k, xResG2, yResG2, zResG2);
                    phi_dx[idxOnG2] = 0.5*(phi[idxXrigt]-phi[idxXleft])/dx;
                }else{
                    phi_dx[idxOnG2] = pfluence[idxOnG2]->getValue()[0];
                }
                phi_dy[idxOnG2] = pfluence[idxOnG2]->getValue()[1];
            }
        }


//******************************************************************************************************////




        for( i = 1; i < xRes+1; i++ ){
            for( j = 1; j < yRes+1; j++ ){
                    idxOnG2  = IDX(i, j, k, xResG2, yResG2, zResG2);
                    idxXleft = IDX(i-1, j, k, xResG2, yResG2, zResG2);
                    idxXrigt = IDX(i+1, j, k, xResG2, yResG2, zResG2);
                    phi_dxdx[idxOnG2] = 0.5*(phi_dx[idxXrigt]-phi_dx[idxXleft])/dx;

                    idxYleft = IDX(i, j-1, k, xResG2, yResG2, zResG2);
                    idxYrigt = IDX(i, j+1, k, xResG2, yResG2, zResG2);
                    phi_dydy[idxOnG2] = 0.5*(phi_dy[idxYrigt]-phi_dy[idxYleft])/dy;
                    phi_dxdy[idxOnG2] = 0.5*(phi_dx[idxYrigt]-phi_dx[idxYleft])/dy;
            }
        }
        for( int ijkG2 = 0; ijkG2 < nG2; ijkG2++ ){
            gridMng->setVectorVariableForNodeG2(ijkG2, AUXVAR2, 0, phi_dxdx[ijkG2]);
            gridMng->setVectorVariableForNodeG2(ijkG2, AUXVAR2, 1, phi_dydy[ijkG2]);
            gridMng->setVectorVariableForNodeG2(ijkG2, AUXVAR2, 2, phi_dxdy[ijkG2]);
         }

        gridMng->sendBoundary2Neighbor(AUXVAR2);

        for( int ijkG2 = 0; ijkG2 < nG2; ijkG2++ ){
            phi_dxdx[ijkG2] = aux2[ijkG2]->getValue()[0];
            phi_dydy[ijkG2] = aux2[ijkG2]->getValue()[1];
            phi_dxdy[ijkG2] = aux2[ijkG2]->getValue()[2];
        }


//****** serve box boundaries ******////
        if (loader->neighbors2Send[NEIGHBOR_LEFT] == MPI_PROC_NULL){

            for( j = 0; j < yResG2; j++ ){
                idxOnG2 = IDX(0, j, k, xResG2, yResG2, zResG2);
                phi_dxdy[idxOnG2] = 0.0;
            }

            idxOnG2 = IDX(0, 0, k, xResG2, yResG2, zResG2);
            idxXrigt = IDX(1, 0, k, xResG2, yResG2, zResG2);
            phi_dxdx[idxOnG2] = (phi_dx[idxXrigt]-phi_dx[idxOnG2])/dx;

            idxYrigt = IDX(0, 1, k, xResG2, yResG2, zResG2);
            phi_dydy[idxOnG2] = (phi_dy[idxYrigt]-phi_dy[idxOnG2])/dy;
            
            idxOnG2  = IDX(0, yResG2-1, k, xResG2, yResG2, zResG2);           
            idxXrigt = IDX(1, yResG2-1, k, xResG2, yResG2, zResG2);
            phi_dxdx[idxOnG2] = (phi_dx[idxXrigt]-phi_dx[idxOnG2])/dx;

            idxYleft = IDX(0, yResG2-2, k, xResG2, yResG2, zResG2);
            phi_dydy[idxOnG2] = (phi_dy[idxOnG2]-phi_dy[idxYleft])/dy;

            for( j = 1; j < yResG2-1; j++ ){
                idxOnG2 = IDX(0, j, k, xResG2, yResG2, zResG2);
                idxYleft = IDX(0, j-1, k, xResG2, yResG2, zResG2);
                idxYrigt = IDX(0, j+1, k, xResG2, yResG2, zResG2);
                phi_dydy[idxOnG2] = 0.5*(phi_dy[idxYrigt]-phi_dy[idxYleft])/dy;

                idxXrigt = IDX(1, j, k, xResG2, yResG2, zResG2);
                phi_dxdx[idxOnG2] = (phi_dx[idxXrigt]-phi_dx[idxOnG2])/dx;                
            }
        }


        if (loader->neighbors2Send[NEIGHBOR_RIGHT] == MPI_PROC_NULL){

            for( j = 0; j < yResG2; j++ ){
                idxOnG2 = IDX(xResG2-1, j, k, xResG2, yResG2, zResG2);
                phi_dxdy[idxOnG2] = 0.0;
            }

            idxOnG2 = IDX(xResG2-1, 0, k, xResG2, yResG2, zResG2);
            idxXleft = IDX(xResG2-2, 0, k, xResG2, yResG2, zResG2);
            phi_dxdx[idxOnG2] = (phi_dx[idxOnG2]-phi_dx[idxXleft])/dx;

            idxYrigt = IDX(xResG2-1, 1, k, xResG2, yResG2, zResG2);
            phi_dydy[idxOnG2] = (phi_dy[idxYrigt]-phi_dy[idxOnG2])/dy;

            idxOnG2  = IDX(xResG2-1, yResG2-1, k, xResG2, yResG2, zResG2);
            idxXleft = IDX(xResG2-2, yResG2-1, k, xResG2, yResG2, zResG2);
            phi_dxdx[idxOnG2] = (phi_dx[idxOnG2]-phi_dx[idxXleft])/dx;

            idxYleft = IDX(xResG2-1, yResG2-2, k, xResG2, yResG2, zResG2);
            phi_dydy[idxOnG2] = (phi_dy[idxOnG2]-phi_dy[idxYleft])/dy;


            for( j = 1; j < yResG2-1; j++ ){
                idxOnG2 = IDX(xResG2-1, j, k, xResG2, yResG2, zResG2);                
                idxYleft = IDX(xResG2-1, j-1, k, xResG2, yResG2, zResG2);
                idxYrigt = IDX(xResG2-1, j+1, k, xResG2, yResG2, zResG2);
                phi_dydy[idxOnG2] = 0.5*(phi_dy[idxYrigt]-phi_dy[idxYleft])/dy;

                idxXleft = IDX(xResG2-2, j, k, xResG2, yResG2, zResG2);
                phi_dxdx[idxOnG2] = (phi_dx[idxOnG2]-phi_dx[idxXleft])/dx;
                
            }
        }


        if (loader->neighbors2Send[NEIGHBOR_BOTTOM] == MPI_PROC_NULL){
            for( i = 0; i < xResG2; i++ ){
                idxOnG2 = IDX(i, 0, k, xResG2, yResG2, zResG2);
                phi_dxdy[idxOnG2] = 0.0;
            }

            idxOnG2 = IDX(0, 0, k, xResG2, yResG2, zResG2);
            idxXrigt = IDX(1, 0, k, xResG2, yResG2, zResG2);
            phi_dxdx[idxOnG2] = (phi_dx[idxXrigt]-phi_dx[idxOnG2])/dx;

            idxYrigt = IDX(0, 1, k, xResG2, yResG2, zResG2);
            phi_dydy[idxOnG2] = (phi_dy[idxYrigt]-phi_dy[idxOnG2])/dy; 

            idxOnG2 = IDX(xResG2-1, 0, k, xResG2, yResG2, zResG2);
            idxXleft = IDX(xResG2-2, 0, k, xResG2, yResG2, zResG2);
            phi_dxdx[idxOnG2] = (phi_dx[idxOnG2]-phi_dx[idxXleft])/dx;

            idxYrigt = IDX(xResG2-1, 1, k, xResG2, yResG2, zResG2);
            phi_dydy[idxOnG2] = (phi_dy[idxYrigt]-phi_dy[idxOnG2])/dy;


            for( i = 1; i < xResG2-1; i++ ){
                idxOnG2  = IDX(i, 0, k, xResG2, yResG2, zResG2);
                idxXleft = IDX(i-1, 0, k, xResG2, yResG2, zResG2);
                idxXrigt = IDX(i+1, 0, k, xResG2, yResG2, zResG2);
                phi_dxdx[idxOnG2] = 0.5*(phi_dx[idxXrigt]-phi_dx[idxXleft])/dx;

                idxYrigt = IDX(i, 1, k, xResG2, yResG2, zResG2);
                phi_dydy[idxOnG2] = (phi_dy[idxYrigt]-phi_dy[idxOnG2])/dy;
            }
        }

        if (loader->neighbors2Send[NEIGHBOR_TOP] == MPI_PROC_NULL){

            for( i = 0; i < xResG2; i++ ){
                idxOnG2 = IDX(i, yResG2-1, k, xResG2, yResG2, zResG2);
                phi_dxdy[idxOnG2] = 0.0;
            }

            idxOnG2 = IDX(0, yResG2-1, k, xResG2, yResG2, zResG2);
            idxXrigt = IDX(1, yResG2-1, k, xResG2, yResG2, zResG2);
            phi_dxdx[idxOnG2] = (phi_dx[idxXrigt]-phi_dx[idxOnG2])/dx;

            idxYleft = IDX(0, yResG2-2, k, xResG2, yResG2, zResG2);
            phi_dydy[idxOnG2] = (phi_dy[idxOnG2]-phi_dy[idxYleft])/dy;

            idxOnG2 = IDX(xResG2-1, yResG2-1, k, xResG2, yResG2, zResG2);
            idxXleft = IDX(xResG2-2, yResG2-1, k, xResG2, yResG2, zResG2);
            phi_dxdx[idxOnG2] = (phi_dx[idxOnG2]-phi_dx[idxXleft])/dx;

            idxYleft = IDX(xResG2-1, yResG2-2, k, xResG2, yResG2, zResG2);
            phi_dydy[idxOnG2] = (phi_dy[idxOnG2]-phi_dy[idxYleft])/dy;



            for( i = 1; i < xResG2-1; i++ ){   
                idxOnG2 = IDX(i, yResG2-1, k, xResG2, yResG2, zResG2);
                idxXleft = IDX(i-1, yResG2-1, k, xResG2, yResG2, zResG2);
                idxXrigt = IDX(i+1, yResG2-1, k, xResG2, yResG2, zResG2);
                phi_dxdx[idxOnG2] = 0.5*(phi_dx[idxXrigt]-phi_dx[idxXleft])/dx;

                idxYleft = IDX(i, yResG2-2, k, xResG2, yResG2, zResG2);
                phi_dydy[idxOnG2] = (phi_dy[idxOnG2]-phi_dy[idxYleft])/dy;
                
            }
        }

//******************************************************************************************************////




        double dt=loader->getTimeStep();
        int xTotSize = loader->totPixelsPerBoxSide[0]+1;
        int yTotSize = loader->totPixelsPerBoxSide[1]+1;
        for( i = 0; i < xResG2 ; i++ ){
            for( j = 0; j < yResG2 ; j++ ){

                idxOnG2 = IDX(i ,j ,k , xResG2, yResG2, zResG2);

                det[idxOnG2] = phi_dxdx[idxOnG2]*phi_dydy[idxOnG2] - pow(phi_dxdy[idxOnG2],2);
                det[idxOnG2] = det[idxOnG2] == 0.0 ? 1.0 : det[idxOnG2];
                double posX = phi_dx[idxOnG2]/dx;
                double posY = phi_dy[idxOnG2]/dy;

                int i_close_glob = int(posX);
                int j_close_glob = int(posY);

                i_close_glob = i_close_glob < 0 ? 0 : i_close_glob;
                j_close_glob = j_close_glob < 0 ? 0 : j_close_glob;
                i_close_glob = i_close_glob >= xTotSize ? (xTotSize-1) : i_close_glob;
                j_close_glob = j_close_glob >= yTotSize ? (yTotSize-1) : j_close_glob;

                double x_diff = posX - int(posX);
                double y_diff = posY - int(posY);

                int idx11 = j_close_glob + xTotSize*i_close_glob;
                int idx21 = j_close_glob + xTotSize*(i_close_glob+1); 
                int idx12 = j_close_glob+1 + xTotSize*i_close_glob;
                int idx22 = j_close_glob+1 + xTotSize*(i_close_glob+1); 

                double flu11 = initialProtonFluence[idx11];
                double flu21 = initialProtonFluence[idx21];
                double flu12 = initialProtonFluence[idx12];
                double flu22 = initialProtonFluence[idx22];

                double interpolatedFlu = 0.0;
                interpolatedFlu += (1.0-y_diff)*( (1.0-x_diff)*flu11 + x_diff*flu21 );
                interpolatedFlu +=      y_diff *( (1.0-x_diff)*flu12 + x_diff*flu22 );

                double logval = interpolatedFlu*abs(det[idxOnG2])/bkgrndProtonFluence[idx11];

                if( logval <= 0.0  ){
                    string msg ="[EquationSolver] logval ="
                        +to_string(logval)+"; interpolatedFlu = "+to_string(interpolatedFlu)
                        +"; abs(det[idxOnG2]) = "+to_string(abs(det[idxOnG2])*1000000000)
                        +"; phi_dxdx[idxOnG2] = "+to_string(phi_dxdx[idxOnG2]*1000000000)
                        +"; phi_dydy[idxOnG2] = "+to_string(phi_dydy[idxOnG2]*1000000000)
                        +"; phi_dxdy[idxOnG2] = "+to_string(phi_dxdy[idxOnG2]*1000000000)
                        +"; i = "+to_string(i)
                        +"; j = "+to_string(j);
                    logger->writeMsg(msg.c_str(), DEBUG);
                    logval = 1.0;
                }
                double damp = dampingCoeff[idxOnG2];
                phi[idxOnG2] += dt*log(logval)*damp;
                gridMng->setVectorVariableForNodeG2(idxOnG2, POTENTIAL, 0, phi[idxOnG2]);
            }
        }

        gridMng->sendBoundary2Neighbor(POTENTIAL);

        for( int ijkG2 = 0; ijkG2 < nG2; ijkG2++ ){
            gridMng->setVectorVariableForNodeG2(ijkG2, PROTONFLUENCE, 3, pot[ijkG2]->getValue()[0]);
            gridMng->setVectorVariableForNodeG2(ijkG2, PROTONFLUENCE, 4, phi_dx[ijkG2]);
            gridMng->setVectorVariableForNodeG2(ijkG2, PROTONFLUENCE, 5, phi_dy[ijkG2]);

            gridMng->setVectorVariableForNodeG2(ijkG2, OUTPUTVAR, 0, phi_dx[ijkG2]-pfluence[ijkG2]->getValue()[0]);
            gridMng->setVectorVariableForNodeG2(ijkG2, OUTPUTVAR, 1, phi_dy[ijkG2]-pfluence[ijkG2]->getValue()[1]);                        
        }





    delete[] phi;
    delete[] phi_dx;
    delete[] phi_dy;
    delete[] phi_dxdx;
    delete[] phi_dxdy;
    delete[] phi_dydy;
    delete[] det;

    }catch(...){
        return SOLVE_FAIL;
    }

    if (rank == 0){
        string msg ="[EquationSolver] SOLVER step ="+to_string(time);
        logger->writeMsg(msg.c_str(), DEBUG);
    }
    
    return SOLVE_OK;
}


EquationSolver::~EquationSolver(){
    finilize();
    logger->writeMsg("[EquationSolver] delete...OK", DEBUG);
}

void EquationSolver::finilize(){
    delete[] initialProtonFluence;
    delete[] dampingCoeff;
}
