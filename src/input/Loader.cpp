#include "Loader.hpp"

using namespace std;

Loader::Loader(){
    logger.writeMsg("[Loader] Loader start!\n", DEBUG);
}

const char  *INIT_CLASS_NAME = "Initializer";
const string  BRACKETS ="()";
const string  BRACKETS_3DOUBLE = "(ddd)";
const string  GET = "get";
const string  GET_TIMESTEP = "getTimestep";
const string  GET_MAX_TIMESTEPS_NUM = "getMaxTimestepsNum";

const string  GET_TIMESTEP_WRITE = "getOutputTimestep";
const string  GET_OUTPUT_DIR = "getOutputDir";

const string  GET_FILENAME_TEMPLATE = "getOutputFilenameTemplate";

const string  LEFT  = "left";
const string  RIGHT  = "right";
const string  GET_DAMPING_BOUNDARY_WIDTH  = "getDampingBoundaryWidth";

const string  GET_FLUENCE = "getProtonFluence";
const string  GET_BKGRD_FLUENCE = "getBackgroundProtonFluence";

const string  GET_MPI_DOMAIN_NUM = "mpiDomainNum";



const string  dirs[] = {"X", "Y", "Z"};

void Loader::load()
{
        MPI_Comm com;
        int periods[3]; //1 - periodic for MPI
        periods[0] = 0;
        periods[1] = 0;
        periods[2] = 1;
    
        int rank ;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        int cs[3];
        int ss[3];


        string msg = "[Loader] initialization...for proc = " + to_string(rank);
        logger.writeMsg(msg.c_str(), DEBUG);
        Py_Initialize();
        this->pInstance = getPythonClassInstance(INIT_CLASS_NAME);


        for (int n=0; n<3; n++){
             string domainNum =GET+dirs[n]+GET_MPI_DOMAIN_NUM;
             this->mpiDomains[n] = (int) callPyLongFunction( pInstance, domainNum, BRACKETS );
        }
    
        MPI_Cart_create(MPI_COMM_WORLD, 3, mpiDomains, periods, 1, &com);
        MPI_Cart_coords(com, rank, 3, cs);
        for(int i = 0; i < 3; i++){
            this->mpiCoords[i] = cs[i];
        }

        int neighborRank;
        for(int i = 0; i < 27; i++){
            this->neighbors2Send.push_back(MPI_PROC_NULL);
            this->neighbors2Recv.push_back(MPI_PROC_NULL);
        }
        for (int a = -1; a <= +1; a++){
            for (int b = -1; b <= +1; b++){
                for (int c = -1; c <= +1; c++){
                        /* __ guess the coordinate of the neighbor __ */
                        ss[0] = cs[0]+a;
                        ss[1] = cs[1]+b;
                        ss[2] = cs[2]+c;

                        if ( (ss[0] < 0 || ss[0] >= mpiDomains[0] ) 
                            || ( ss[1] < 0 || ss[1] >= mpiDomains[1] ) 
                                || ( ss[2] < 0 || ss[2] >= mpiDomains[2] ) ){
                        
                            neighborRank = MPI_PROC_NULL;
                        }else{
                            MPI_Cart_rank(com, ss, &neighborRank);
                        }
                        this->neighbors2Send[(1+c)+3*((1+b)+3*(1+a))] = neighborRank;
                        this->neighbors2Recv[(1-c)+3*((1-b)+3*(1-a))] = neighborRank;
                    }
                }
            }


        double boxSizePerDomain[3];
        for (int n=0; n<3; n++){
        
            string res = GET + dirs[n] + "resolution";
            string right_str= GET + dirs[n] + "right";
            this->boxCoordinates[n][0] = 0.0; //origin point (0,0,0) is left lower coner i=0 j=0 k=0
            this->boxCoordinates[n][1] = callPyFloatFunction( pInstance, right_str, BRACKETS );
            
            
            this->boxSizes[n] = boxCoordinates[n][1];
            //  length of the box in normalized units
            
            //  length of the domain in normalized units
            boxSizePerDomain[n] = boxSizes[n]/mpiDomains[n];
            
            // this->boxCoordinates[n][0] = boxCoordinates[n][0] + cs[n]*boxSizePerDomain[n];
            // this->boxCoordinates[n][1] = boxCoordinates[n][0] + boxSizePerDomain[n];

            //  number of pixel per box
            this->totPixelsPerBoxSide[n] = (int) callPyLongFunction( pInstance, res, BRACKETS);
            
            if( boxSizes[n] <= 0.0 ){
                throw runtime_error("box has zero length!");
            }

           //  number of pixel per box
            this->totPixelsPerBoxSide[n] = (int) callPyLongFunction( pInstance, res, BRACKETS);
            
            if( totPixelsPerBoxSide[n] == 1 ){
                // length per pixel equals to total size
                this->spatialSteps[n]=boxSizes[n];
            }else{
                this->spatialSteps[n]=boxSizes[n]/(double)(totPixelsPerBoxSide[n]);
            }
            
            double approxRes = int(totPixelsPerBoxSide[n]/double(mpiDomains[n]));
            
            this->resolution[n] = int(approxRes);
            
            int delta = totPixelsPerBoxSide[n] - resolution[n]*mpiDomains[n];
            
            this->offsetInPixels[n] = 0;
            
            if ( delta == 0.0 ){
                //  length of the domain in normalized units
                boxSizePerDomain[n]     = resolution[n]*spatialSteps[n];
                this->offsetInPixels[n] = resolution[n]*cs[n];// global offset in pixels
                
            } else if ( cs[n] < delta ){
                this->resolution[n] += 1;
                
                boxSizePerDomain[n]     = resolution[n]*spatialSteps[n];
                this->offsetInPixels[n] = resolution[n]*cs[n];
                
            } else {
                boxSizePerDomain[n]     = resolution[n]*spatialSteps[n];
                this->offsetInPixels[n] = (resolution[n]+1)*delta + resolution[n]*(cs[n]-delta);
            }
            
            this->boxCoordinates[n][0] = offsetInPixels[n]*spatialSteps[n];
            this->boxCoordinates[n][1] = boxCoordinates[n][0]+boxSizePerDomain[n];

    }


    for (int n=0; n<2; n++){
            string right_damping = GET_DAMPING_BOUNDARY_WIDTH + dirs[n] + RIGHT;
            string left_damping  = GET_DAMPING_BOUNDARY_WIDTH + dirs[n] + LEFT;
            
            PyObject* calldampingBoundaryWidthMethod = getPyMethod( pInstance, left_damping, BRACKETS );
            
            if( calldampingBoundaryWidthMethod != NULL ){
                this->dampingBoundaryWidth[n][0] = PyFloat_AsDouble(calldampingBoundaryWidthMethod);
            }
            
            calldampingBoundaryWidthMethod = getPyMethod( pInstance, right_damping, BRACKETS );
            
            if( calldampingBoundaryWidthMethod != NULL ){
                this->dampingBoundaryWidth[n][1] = PyFloat_AsDouble(calldampingBoundaryWidthMethod);
            }        
    }
    
    this->timeStep           = callPyFloatFunction( pInstance, GET_TIMESTEP, BRACKETS );
        
    this->maxTimestepsNum    = callPyFloatFunction( pInstance, GET_MAX_TIMESTEPS_NUM, BRACKETS );
    
    this->timestepsNum2Write = callPyFloatFunction( pInstance, GET_TIMESTEP_WRITE, BRACKETS );
    
    this->outputDir          = callPyStringFunction( pInstance, GET_OUTPUT_DIR, BRACKETS );
    
    this->fileNameTemplate   = callPyStringFunction( pInstance, GET_FILENAME_TEMPLATE, BRACKETS );
        
    
    if(rank == 0){


        printf("[Loader]  Box size: \n");

        for (int n=0; n<3; n++){
            printf("[Loader]  [%s] [%1.5f, %1.5f] l = %1.5f res = %3i => step  = %1.6f\n",
                   dirs[n].c_str(), boxCoordinates[n][0],boxCoordinates[n][1],
                   boxSizes[n],resolution[n],spatialSteps[n]);
        }
        
        
        msg = "[Loader] timeStep = "+to_string(timeStep);
        logger.writeMsg(msg.c_str(), INFO);
        msg = "[Loader] Max timeStep number = "+to_string(maxTimestepsNum);
        logger.writeMsg(msg.c_str(), INFO);
        msg = "[Loader]  timeStep number to write to file = "+to_string(timestepsNum2Write);
        logger.writeMsg(msg.c_str(), INFO);
        msg = "[Loader]  outputDir = "+outputDir;
        logger.writeMsg(msg.c_str(), INFO);
        msg = "[Loader]  filename template = "+fileNameTemplate;
        logger.writeMsg(msg.c_str(), INFO);

        if( this->dampingBoundaryWidth[0][0]>0 || this->dampingBoundaryWidth[0][1]>0 ){
            string dampingBCmsg = "[Loader] [BC - X] damping layer width: left = "
                                +to_string( this->dampingBoundaryWidth[0][0] )
                                +" / right = "+to_string(this->dampingBoundaryWidth[0][1] );
            logger.writeMsg(dampingBCmsg.c_str(), INFO);
        }
    
        if( this->dampingBoundaryWidth[1][0]>0 || this->dampingBoundaryWidth[1][1]>0 ){
            string dampingBCmsg = "[Loader] [BC - Y] damping layer width: left = "
                                +to_string( this->dampingBoundaryWidth[1][0] )
                                +" / right = "+to_string(this->dampingBoundaryWidth[1][1] );
            logger.writeMsg(dampingBCmsg.c_str(), INFO);
        }

        logger.writeMsg("[Loader] initialization...OK!", DEBUG);
    }

}

double Loader::getTimeStep(){
    return timeStep;
}

int Loader::getMaxTimestepsNum(){
    return maxTimestepsNum;
}

int Loader::getTimestepsNum2Write(){
    return timestepsNum2Write;
}

string Loader::getOutputDir() const {
    return outputDir;
}

string Loader::getFilenameTemplate(){
    return fileNameTemplate;
}

double Loader::getProtonFluence(double x, double y, double z){
    return callPyFloatFunctionWith3args( pInstance, GET_FLUENCE, BRACKETS_3DOUBLE, x, y, z );
}

double Loader::getBackgroundProtonFluence(double x, double y, double z){
    return callPyFloatFunctionWith3args( pInstance, GET_BKGRD_FLUENCE, BRACKETS_3DOUBLE, x, y, z );
}


Loader::~Loader(){
    Py_Finalize();
    logger.writeMsg("[Loader] FINALIZE...OK!", DEBUG);
}


/*#################################################################################################*/


PyObject * Loader::getPythonClassInstance(string className){
    PyObject  *pName, *pClass, *pModule, *pDict;
    string msg = "[Loader] Start to instantiate the Python class " + className;
    logger.writeMsg(msg.c_str(), DEBUG);
    
    #if PY_MAJOR_VERSION >= 3
    pName = PyUnicode_FromString(className.c_str());
    #else
    pName = PyString_FromString(className.c_str());
    #endif
    
    pModule = PyImport_Import(pName);
    
    if( pModule == NULL ){
        logger.writeMsg("*****************************************************", CRITICAL);
        logger.writeMsg("****                                             ****", CRITICAL);
        logger.writeMsg("****  STOP SIMULATION!!!    INPUT FILE PROBLEM   ****", CRITICAL);
        logger.writeMsg("****                                             ****", CRITICAL);
        logger.writeMsg("****   try debug mode 'make -DLOG'               ****", CRITICAL);
        logger.writeMsg("****   check indents in python input file        ****", CRITICAL);
        logger.writeMsg("****   check python enviroment and imports       ****", CRITICAL);
        logger.writeMsg("*****************************************************", CRITICAL);
        exit(-1);
    }
    
    pDict = PyModule_GetDict(pModule);
    pClass = PyDict_GetItemString(pDict, className.c_str());
    
    if( PyCallable_Check(pClass) ){
        pInstance = PyObject_CallObject(pClass, NULL);
    }else{
        logger.writeMsg("[Loader] Cannot instantiate the Python class", CRITICAL);
        pInstance = nullptr;
    }
    
    logger.writeMsg("[Loader] finish to instantiate the Python class ", DEBUG);
    
    return pInstance;
}
    


double Loader::callPyFloatFunction( PyObject* instance,
                                   const string funcName,
                                   const string brackets){
    return PyFloat_AsDouble(getPyMethod(instance,funcName,brackets));
}

double Loader::callPyFloatFunctionWith3args( PyObject* instance,
                                            const string funcName,
                                            const string brackets,
                                            double x,double y,double z){
    
    PyObject *methodReturn = PyObject_CallMethod(instance, strdup(funcName.c_str()),
                                 strdup(BRACKETS_3DOUBLE.c_str()),x,y,z);
    
    if( methodReturn ){
        return PyFloat_AsDouble(methodReturn);
    }else{
        PyErr_Print();
        throw runtime_error("Problem to call python function: "+funcName);
    }
    
}

long Loader::callPyLongFunction( PyObject* instance,
                                const string funcName,
                                const string brackets){
    
    return PyLong_AsLong(getPyMethod(instance,funcName,brackets));
}

string Loader::callPyStringFunction( PyObject* instance,
                                    const string funcName,
                                    const string brackets){
    
    #if PY_MAJOR_VERSION >= 3
    PyObject* str = PyUnicode_AsUTF8String(getPyMethod(instance,funcName,brackets));
    return PyBytes_AsString(str);
    #else
    return PyString_AS_STRING(getPyMethod(instance,funcName,brackets));
    #endif

}

PyObject* Loader::getPyMethod(PyObject* instance,
                              const string funcName,
                              const string brackets){
    
    return PyObject_CallMethod(instance, strdup(funcName.c_str()),
                               strdup(brackets.c_str()));
}

/*#################################################################################################*/
