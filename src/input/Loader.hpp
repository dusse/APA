#ifndef Loader_hpp
#define Loader_hpp

#include <Python.h>
#include <mpi.h>
#include <stdio.h>
//#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <math.h>
#include "../misc/Logger.hpp"
#include "../misc/Misc.hpp"

class Loader{
    
private:
    Logger logger;
    double timeStep;
    
    int maxTimestepsNum;
    int timestepsNum2Write;
    
    std::string outputDir;
    std::string fileNameTemplate;

    PyObject *pInstance;
    
    double callPyFloatFunction( PyObject*, const std::string, const std::string);
    
    double callPyFloatFunctionWith3args( PyObject*, const std::string, const std::string,
                                                double, double, double);
    
    long callPyLongFunction( PyObject*, const std::string, const std::string);
    
    std::string callPyStringFunction( PyObject*, const std::string, const std::string);
    
    PyObject* getPyMethod(PyObject* , const std::string , const std::string);
    
public:
    
    int resolution[3];
    int totPixelsPerBoxSide[3];
    int offsetInPixels[3];
    int mpiDomains[3];
    int mpiCoords[3];

    double boxCoordinates[3][2];
    double boxSizes[3];
    double spatialSteps[3];

    double dampingBoundaryWidth[2][2];
    
       //MPI staff
    std::vector<int> neighbors2Send;//27
    std::vector<int> neighbors2Recv;//27
    
    
    Loader();
    ~Loader();
    void load();
    double getProtonFluence(double,double,double);
    double getBackgroundProtonFluence(double,double,double);
    double getTimeStep();
    int getMaxTimestepsNum();
    int getTimestepsNum2Write();
    std::string getOutputDir() const;
    std::string getFilenameTemplate();
    PyObject * getPythonClassInstance(std::string className);
    
};
#endif /* Loader_hpp */
