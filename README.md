# QonGPU

Program from my Bachelorthesis. Its main purpose is to simulate a time resolved 1D ionization with some
python tools to analyze the data.

## Method
For any given arbitrary potential

##Dependencies
1. CMake >= 3.0
2. Boost >= 1.5.6
3. GCC >= 4.9.2
4. HDF5/HDF5lite >= 1.86
5. Python 3: Matplotlib & Numpy for data analysis
6. CUDA comptible with GCC >= 4.9.2

##Building
Since this is single-node HPC application only Linux is currently supported. 

```{r, engine='bash', count_lines}
cd /location_to_QonGPU/
./configure.sh

```
The script will take care of all the build folders etc. The results are saved per default
in a HDF5 file called  ```  res.h5``` and is located in the build folder. The way the timesteps are saved is that
one timesteps n creates  2 HDF5 datasets, one for the imaginary and one for the real part of the
wave-function. These datasets will always have the name ``` dsetNreal``` and ```dsetNimg```, where
N is the number of timesteps calculated. The stationary states as well as the energies of each stationary state are saved inside the
HDF5 file ```build/sim1.h5```. All the simulation parameters are set inside the ```main.cu``` file.
