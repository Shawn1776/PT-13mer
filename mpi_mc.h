#include "mpi.h"

using namespace MPI; 
class mc_mpi{
public:
	void mpi_init();



};

void mpi_mc::mpi_init(int &totalid, int &id){
	//MPI::Init ( argc, argv );
	MPI::Init();
	totalid = MPI::COMM_WORLD.Get_size();
	id = MPI::COMM_WORLD.Get_rank();
}