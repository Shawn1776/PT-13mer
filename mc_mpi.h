#include "mpi.h"

//using namespace MPI; 
class mc_mpi{
public:
	void repilca_exchange(int iter_group_id,int myid,MPI_COMM_WORLD);



};

	void mc_mpi::energy_exchange(int iter_group_id,int myid,MPI_COMM_WORLD, double E, int &local_flag){  // post - condition --> change the local_flag
// this function has 2 parts, part1-> exchange then, decide if swap repila
		// exchange
		double localEnergy = E, receivedEnergy = 0.0;
		int currentid = id, otherid = -1;
		// int local_flag=0, recv_flag=0;
		
		if ( iter_group_id%2==0 ){ // e.g. iter_group_id = 0, 2, 4 ~~~
			if( myid%2==0 ){ // myid = 0, 2, 4 ~~~
			    MPI::COMM_WORLD.Sendrecv(/*from myid send to myid+1*/&localEnergy,1, MPI::DOUBLE, (myid+1),/*tag1*/myid,
				   		 		     /*form myid+1 send to myid*/,&receivedEnergy,1, MPI::DOUBLE, (myid+1),/*tag2*/(myid+1), 
		                		status);
				// PT crertia
				double p_accept = std::min( 1.0, exp((localEnergy-receivedEnergy)*(1.0/Ts[myid] - 1.0/Ts[myid+1])) );
				double r = RAN01();
				// send the flag
				if ( p_accept>r ){
					// repicla_exchange
					local_flag = 1;
					MPI::COMM_WORLD.Send(&local_flag, 1, MPI::INT, (myid+1),myid*500);
				}else{
					// send the flag = 0
					MPI::COMM_WORLD.Send(&local_flag, 1, MPI::INT, (myid+1),myid*500);
					cout<<"myid: "<<myid<<";  local_flag is: "<<local_flag<<endl;
				}//if- else
				if(local_flag){
					MPI::COMM_WORLD.Sendrecv(/*from myid send to myid+1*/&localEnergy,1, MPI::DOUBLE, (myid+1),/*tag1*/myid,
				   		 		     /*form myid+1 send to myid*/,&receivedEnergy,1, MPI::DOUBLE, (myid+1),/*tag2*/(myid+1), 
		                		status);
				}

			}else{			 // myid = 1, 3, 5 ~~~
				Comm::Sendrecv(/*from myid send to myid-1*/&localEnergy,1, MPI::DOUBLE, (myid-1),/*tag3*/myid,
				 		   /*form myid-1 send to myid*/,&receivedEnergy,1, MPI::DOUBLE, (myid-1),/*tag4*/myid-1, 
                			status);
					// recv the flag
					MPI::COMM_WORLD.Recv(&local_flag, 1, MPI::INT, (myid-1),(myid-1)*500);

					// 
					cout<<"myid: "<<myid<<";  local_flag is: "<<local_flag<<endl;
			}// if - myid%2---> repilca exchange or not
			MPI::COMM_WORLD.Barrier();

		}

		if (iter_group_id%2!=0){

		}

	}