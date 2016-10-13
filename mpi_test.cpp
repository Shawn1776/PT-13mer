#include "mpi.h"
#include <stdio.h>
#include <vector>
#include "ranMARS.h"
#include <cmath>

using namespace std;
int main(int argc, char *argv[]){
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //    MPI Initialization //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  MPI::Init();
  int totalid = MPI::COMM_WORLD.Get_size();
  int myid = MPI::COMM_WORLD.Get_rank();
  init_RAN(22+myid,34,56,78); 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //    Temperature Initialization //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  vector<double> Ts{1.000000,1.071429,1.153846, 1.250000};
  double localEnergy, receivedEnergy;
  int local_flag = -1;
  MPI::Status status;

    if( myid%2==0 ){ 
      localEnergy=-1.0; 
      receivedEnergy=0.0;
     // cout<<"Before--myid: "<<myid<<";  localEnergy is: "<<localEnergy<<" recvEnergy: "<<receivedEnergy<<endl;
    }else{
      localEnergy=1.0; 
      receivedEnergy=0.0;
      // cout<<"Before--myid: "<<myid<<";  localEnergy is: "<<localEnergy<<" recvEnergy: "<<receivedEnergy<<endl;
    }

    // int local_flag=0, recv_flag=0;
    
    // if ( iter_group_id%2==0 ){ // e.g. iter_group_id = 0, 2, 4 ~~~
      if( myid%2==0 ){ // myid = 0, 2, 4 ~~~
          cout<<"before--myid: "<<myid<<";  local_flag is: "<<local_flag<<endl;
          MPI::COMM_WORLD.Sendrecv(/*from myid send to myid+1*/&localEnergy,1, MPI::DOUBLE, (myid+1),/*tag1*/myid,
                       /*form myid+1 send to myid*/&receivedEnergy,1, MPI::DOUBLE, (myid+1),/*tag2*/(myid+1), 
                        status);
        // PT crertia
        // double p_accept = std::min( 1.0, exp((localEnergy-receivedEnergy)*(1.0/Ts[myid] - 1.0/Ts[myid+1])) );
        double p_accept = 1.0;
        double r = RAN01();
        // send the flag
        if ( p_accept>r ){
          // repicla_exchange
          local_flag = 1;
          MPI::COMM_WORLD.Send(&local_flag, 1, MPI::INT, (myid+1),myid*500);
          cout<<"After--myid: "<<myid<<";  local_flag is: "<<local_flag<<endl;
        }else{
          // send the flag = 0
          local_flag = 0;
          MPI::COMM_WORLD.Send(&local_flag, 1, MPI::INT, (myid+1),myid*500);
          cout<<"After--myid: "<<myid<<";  local_flag is: "<<local_flag<<endl;
        }//if- else
        if(local_flag){
            // exchange repilca 
        }

      }else{       // myid = 1, 3, 5 ~~~
        cout<<"before--myid: "<<myid<<";  local_flag is: "<<local_flag<<endl;
        MPI::COMM_WORLD.Sendrecv(/*from myid send to myid-1*/&localEnergy,1, MPI::DOUBLE, (myid-1),/*tag3*/myid,
               /*form myid-1 send to myid*/&receivedEnergy,1, MPI::DOUBLE, (myid-1),/*tag4*/myid-1, 
                        status);
          // recv the flag
          MPI::COMM_WORLD.Recv(&local_flag, 1, MPI::INT, (myid-1),(myid-1)*500);

          // 
          cout<<"After--myid: "<<myid<<";  local_flag is: "<<local_flag<<endl;
      }// if - myid%2---> repilca exchange or not
      MPI::COMM_WORLD.Barrier();





  MPI::Finalize();


  return 0;
}








