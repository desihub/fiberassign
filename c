if [ "$1" == "d" ]; then
	mv Assign.* outputs
	cp -R src/* ../tmp/src 
	cd src
	make all
	cd ..
	qsub run
	qs -u lgarrigu

elif [ "$1" == "cc" ]; then
	cd src
	#make clean
	make all

elif [ "$1" == "c" ]; then
	mv Assign.* outputs
	cp -R src/* ../tmp/src 

elif [ "$1" == "v" ]; then
	valgrind ./assign /project/projectdirs/desi/mocks/preliminary/objects0.rdzipn /project/projectdirs/desi/software/edison/desimodel/0.3.1/data/footprint/desi-tiles.par /project/projectdirs/desi/software/edison/desimodel/0.3.1/data/focalplane/fiberpos.txt assignment

elif [ "$1" == "q" ]; then
	./assign /project/projectdirs/desi/mocks/preliminary/objects0.rdzipn /project/projectdirs/desi/software/edison/desimodel/0.3.1/data/footprint/desi-tiles.par /project/projectdirs/desi/software/edison/desimodel/0.3.1/data/focalplane/fiberpos.txt assignment

elif [ "$1" == "s" ]; then
	module swap PrgEnv-intel PrgEnv-gnu

fi
