#include "mc_mpi.h"
#include "polymer.h"
#include "ranMARS.h"
#include "MonteCarlo.h"
#include <fstream>

int main(int argc, char *argv[]){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    // MPI Initialization //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	MPI::Init();
//	int totalid = MPI::COMM_WORLD.Get_size();
//	int id = MPI::COMM_WORLD.Get_rank();	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    // Temperature Initialization //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	vector<double> Ts{1.000000,1.071429,1.153846, 1.250000};
//	double temperature= Ts[id];
	double temperature = 1.0;
	int id = 11111;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    // RNG() Initialization //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	init_RAN(22+id,34,56,78);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    // Polymer and MC system Initialization //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	ploymer p1(13);
	MonteCarlo mc1;
	int size = floor((p1.coordv.size()+1)/3);
	cout<<"size: "<<size<<endl;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    // Output files setups //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
    string filename =to_string(id).append("_E_").append(to_string(temperature));
	//filename.append(to_string(temperature));
	ofstream outfile;
	outfile.open(filename, ios::out);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    // All ploymers need to reach equilibrium first //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int equilibration = size*5000;
		cout<<"I am rank "<<id<<endl;


	//for (int t=0; t<equilibration; t++){
		//cout<<"I am rankkkkkkkkkkkk "<<id<<endl;

		double E = mc1.Metropolis_single_monomer(size,p1.coordv,temperature); // efficince, need to pass the calculated E back to the function
			// if ((t+1)%1==0){
			// 	//cout << (t+1)/1300 <<"	"<< E <<endl; // output every 1300, 100 swap
			// 	outfile <<(t+1) <<"	"<< E <<endl;
			// }//if
		//cout<<"I am rankkkkkkkkkkkk "<<id<<endl;	
	//	}//for
cout<<"I am rankkkkk "<<id<<endl;
		outfile.close();
//	MPI::Finalize();
	return 0;
}


