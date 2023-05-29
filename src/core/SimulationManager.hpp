//
//  SimulationManager.hpp

#ifndef SimulationManager_hpp
#define SimulationManager_hpp

#include <stdio.h>
#include <chrono>
#include <mpi.h>
#include "../input/Loader.hpp"
#include "../grid/GridManager.hpp"

#include "../misc/Logger.hpp"
#include "../output/Writer.hpp"

#include "../solvers/EquationSolver.hpp"
#include "../solvers/ModelInitializer.hpp"

#include "../common/variables/VectorVar.hpp"

class SimulationManager {
private:
    std::unique_ptr<Logger> logger;
    std::shared_ptr<GridManager> gridMng;
    
    std::shared_ptr<ModelInitializer> initMng;
        
    std::shared_ptr<EquationSolver> solver;
    
    std::shared_ptr<Loader> loader;
    std::unique_ptr<Writer> writer;

public:
    SimulationManager(int ac, char **av);
    void initialize();
    void runSimulation(int ac, char **av);
    void finilize();
    ~SimulationManager();
};

#endif
