#include "polymer.h"
#include "ranMARS.h"
#include "MonteCarlo.h"
#include <fstream>


using namespace std;
typedef vector<double> vd; // shorthand

int main(){
	int id = 10;
	init_RAN(22+id,34,56,78);
	// intialize
	ploymer p1(13);
	MonteCarlo mc1;

	// must calculate the size to optimize the code
	int size = floor((p1.coordv.size()+1)/3);
	double temperature = 1.0;
	// Mento Carlo Updates, each mono-mer
	int totalMC = 13*100000000; //E5 *000; // 100 MC swap
	string filename ="E_";
	for(int tn =0; tn<4; tn++){ // tn-> temperature number
		double temperature = 1.0;
		temperature += tn*0.1;
		//string filename= static_cast<ostringstream*>( &(ostringstream() << temperature) )->str();
		string filename ="E_";
		filename.append(to_string(temperature));
		ofstream outfile;
		outfile.open(filename, ios::out);
		
		for (int t =0; t<totalMC; t++){
			double E = mc1.Metropolis_single_monomer(size,p1.coordv,temperature); // efficince, need to pass the calculated E back to the function
			if ((t+1)%1300==0){
				//cout << (t+1)/1300 <<"	"<< E <<endl; // output every 1300, 100 swap
				outfile <<(t+1)/1300 <<"	"<< E <<endl;
			}
			
		}
		outfile.close();

	}
	
	//p1.outputCoordv();
	return 0;
}
