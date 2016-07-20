# QonGPU
A project to simulate arbitrary quantum systems on GPUs. Currently the simulation is only able
to hadle one-dimensional systems, but a broadening of spatial degrees of freedom will be done in the future.

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
cd /location_to_QonGPU/
./configure.sh

```
The script will take care of all the build folders etc. The results are saved per default
in a HDF5 file called  ```  res.h5``` and is located in the build folder. The way the timesteps are saved is that
one timesteps n creates  2 HDF5 datasets, one for the imaginary and one for the real part of the
wave-function. These datasets will always have the name ``` dsetN+1real``` and ```dsetN+1img```, where
N is the number of timesteps calculated. The reason for ```N+1``` is to avoid an overwrite of  the first
dataset, this will be changed in the future. The stationary states as well as the energies of each stationary state are saved inside the
HDF5 file ```build/sim1.h5```. All the simulation parameters are set inside the ```main.cu``` file.