#include "Writer.hpp"

using namespace std;
using namespace chrono;

Writer::Writer(shared_ptr<Loader> load, shared_ptr<GridManager> gridMnr):
                loader(move(load)), gridMgr(move(gridMnr)){
    logger.reset(new Logger());

    this->outputDir       = loader->getOutputDir();
    this->fileNamePattern = loader->getFilenameTemplate();
    
    logger->writeMsg("[Writer] initialize...OK", DEBUG);
}

void Writer::write(int fileNum){
    auto start_time = high_resolution_clock::now();
   
    logger->writeMsg("[Writer] writing...", DEBUG);
    int subdomainXSize = loader->resolution[0];
    int subdomainYSize = loader->resolution[1];
    int subdomainZSize = loader->resolution[2];
    
    int size[3] = {subdomainXSize, subdomainYSize, subdomainZSize};
    
    int totalNodeNum = subdomainXSize*subdomainYSize*subdomainZSize;
    int ijNode;
    vector<vector<VectorVar>> vectorVars = gridMgr->getVectorVariablesForAllNodes();
        
    MPI_Info info = MPI_INFO_NULL;
    
    hid_t access = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(access, MPI_COMM_WORLD, info);
    
    string fileName = outputDir + fileNamePattern + to_string(fileNum) + ".h5";
    hid_t fileID = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, access);
    
    const string groupname = "/vars";
    hid_t group   = H5Gcreate(fileID, groupname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
    
    int vars4EachNode = vectorVars[0].size();
    for (int idx_var = 0; idx_var < vars4EachNode; idx_var++) {
        
        int varDim = vectorVars[0][idx_var].getSize();
        for (int dir = 0; dir < varDim; dir++) {
            
        string varName  = to_string(idx_var) +"_"+ to_string(dir);
            double* var = new double[totalNodeNum];
            
            for (ijNode = 0; ijNode < totalNodeNum; ijNode++) {
                var[ijNode] = (vectorVars[ijNode][idx_var].getValue())[dir];
            }
            
            writeParallel(fileID, group, dxpl_id, varName, var, size);
            
            delete [] var;
        }
    }

    
    H5Gclose(group);
    H5Fflush(fileID, H5F_SCOPE_GLOBAL);
    H5Pclose(dxpl_id);
    H5Fclose(fileID);
        
    auto end_time = high_resolution_clock::now();
    string msg ="[Writer] writing duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}




void Writer::writeParallel(hid_t fileID, hid_t group, hid_t dxpl_id,
                           string dsetName , const double* data, int sizes[3]){

    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    hsize_t offset[3];
    hsize_t nLoc[3];
    hsize_t nGlob[3];

    for (int dir = 0; dir < 3; dir++) {
        offset[dir] = loader->offsetInPixels[dir];
        nLoc[dir]   = sizes[dir];
        nGlob[dir]  = loader->totPixelsPerBoxSide[dir];
    }

    hid_t memspace  = H5Screate_simple(3, nLoc , NULL);
    hid_t filespace = H5Screate_simple(3, nGlob, NULL);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, nLoc, NULL);

    hid_t dset;
    dset  = H5Dcreate(group, dsetName.c_str(), H5T_NATIVE_DOUBLE, filespace,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, memspace, filespace, dxpl_id, data);
    H5Dclose(dset);

    H5Sclose(memspace);
    H5Sclose(filespace);
}






