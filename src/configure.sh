
#!/bin/sh


cmake .

p


#check number of cores for later use with make
cores=1
if [ $(which nproc) ]
then
	cores=$(nproc)
elif [ $(which grep) -a -f /proc/cpuinfo ]
then
	cores=$(grep -c ^processor /proc/cpuinfo)
fi

echo "Found number of Cores: ${cores}. Use $(expr ${cores} + 1) jobs for parallel build using make."

echo "Running make"

make -j$((cores+1)) --no-print-directory

echo "Finished"

echo "Running qsolve"
#build/Solver/qsolve



