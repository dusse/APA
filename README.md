
***************************************                                   
         A       PPPPPP       A        
        A A      P    P      A A       
       A   A     P    P     A   A      
      AAAAAAA    PPPPPP    AAAAAAA     
     A       A   P        A       A    
    A         A  P       A         A   
 Archie's Proton-radiography Algorithm 
****************************************

# APA
 Archie's Proton-radiography Algorithm 
 for extracting path-integrated magnetic fields.

 Bott, A. F. A., Graziani, C., Tzeferacos, P., White, T. G., 
 Lamb, D. Q., Gregori, G., & Schekochihin, A. A. (2017). 
 Proton imaging of stochastic magnetic fields. 
 Journal of Plasma Physics, 83(6), 905830614.

 stack: C++11, MPI, HDF5, python

_______________________
#       HOWTO
_______________________
1. before 'make' need to set in makefile
    HDF5_PATH= path to hdf5 lib (last well used 1.10.5)
    MPI_PATH= path to mpi lib (last well used openmpi 9.0.0)
    PYTHON27_INC= path to python include
    PYTHON27_LIB= path to python lib

2. for running default example from src/input/Initializer.py
    mpirun -n 2 apa.exe

3. normally need to set input file path containing Initializer.py
    mpirun -n 2 apa.exe PATH/TO/PYTHON/INPUT/FILE

4. before running need to create output folder and set in Initializer.py

5. for visualization use python notebook in folder NOTEBOOK/



