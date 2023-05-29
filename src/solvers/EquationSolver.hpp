//
//  EquationSolver.hpp

#ifndef EquationSolver_hpp
#define EquationSolver_hpp

#include <chrono>
#include <stdio.h>
#include <iostream>
#include <string>
#include <memory>

#include "../grid/GridManager.hpp"
#include "../input/Loader.hpp"
#include "../misc/Misc.hpp"



const static int  SOLVE_OK   = 0;
const static int  SOLVE_FAIL = 1;

class EquationSolver{
    
private:
    
    std::unique_ptr<Logger> logger;
    
    std::shared_ptr<Loader> loader;
    std::shared_ptr<GridManager> gridMng;

    double* initialProtonFluence;
    double* bkgrndProtonFluence;
    double* dampingCoeff;

    
public:
    EquationSolver(std::shared_ptr<Loader>,
                  std::shared_ptr<GridManager>);
        
    void initialize();
    void initDampingCoeff();
    int solve(double);
    void finilize();
    ~EquationSolver();
};
#endif
