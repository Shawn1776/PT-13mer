mpi_main: mpi_main.cpp
	mpic++ -std=c++11 ranMARS.cpp mpi_main.cpp -o mpi_main
#run:
#	mpirun -np 4 ./mpi
clean:
	#find . -type f | xargs -n 5 touch # fix "make: Warning: File `makefile' has modification time 48 s in the future"
	#rm -rf $(OBJS)
	rm ./mpi_main

