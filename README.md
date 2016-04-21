# SchroedingerSolver
This project aims to solve the Schrödinger equation for arbitraty potentials in up to 3 spacial dimensions on CUDA. 

The current focus is to write very fast code in order to solve more complex systems!
## Method
### 1st spacial dimension:
First a static solver is implented using a model potential. This static solver
uses the numerov algorithm to calculate the staionary states.


##Dependencies
1. CMake >= 3.0
2. Boost >= 1.5.6
3. GCC >= 4.9.2
4. HDF5/HDF5lite >= 1.86
5. Python 3: Matplotlib & Numpy for data analysis

##Building
Only Linux support planned!
```{r, engine='bash', count_lines}
cd /location_to_SchroedingerSolver/
./configure.sh

```
Configure will compile the program and run it for you. The
output file configuration is a bit messy and will be fixed later on.
If you want to see the static results:

To get the Python script going just
```{r, engine='bash', count_lines}
 cp src/pyscripts/Level2.py build
 cd build/
 python Level2.py
```
This is still a very unclean way and will be fixed in the future!

