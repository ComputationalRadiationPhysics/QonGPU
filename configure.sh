
#!/bin/sh

mkdir -p build/

cd build
cmake .. || exit 1




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

make -j$((cores+1)) --no-print-directory || exit 1
make test || exit 1

echo "Finished"

echo "Running qsolve"
loc=$(find -name qsolve)

mv ${loc}  ./build/
build/qsolve

#if [ "$1"=="-py" ]; then
#        echo "Use Python graphics script"
#        python pyscripts/Graphicscript.py
#fi
