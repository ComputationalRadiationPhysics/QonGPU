#!/bin/sh
programName="qsolve"

echo "$(tput setab 7)                         $(tput sgr0)"
echo "$(tput setab 7)  $(tput sgr0) $(tput bold)$(tput setaf 4)build $programName$(tput sgr0) $(tput setab 7)  $(tput sgr0)"
echo "$(tput setab 7)                         $(tput sgr0)"

tput setaf 3

#check for cmake, make
if [ ! $(which cmake) -o ! $(which make) ]
then
	echo "Either cmake or make could not be found."
	exit 4
fi
#check current pwd - try to get to parallelParts/ (first one)
currentPwd="$(pwd)"
while [ ! -e "CMakeLists.txt" -a "$(pwd)" != "/" ]
do
	cd ..
done
if [ "$(pwd)" = "/" ]
then
	echo "Please go into ${programName}-directory or in a sub-directory of it."
	exit 3
fi

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

addOption=""

echo "Running cmake"
tput sgr0

cmake "${addOption}" .

tput setaf 3
echo "Running make"
tput sgr0

make -j$((cores+1)) --no-print-directory

tput setaf 3
echo "Finished"
tput sgr0
