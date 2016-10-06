#include "polymer.h"
#include "ranMARS.h"
#include "MonteCarlo.h"


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
	int totalMC = 13*1000000; // 100 MC swap
	for (int t =0; t<totalMC; t++){
		if ((t+1)%13==0){
			cout << (t+1)/13 <<"	"<<mc1.Metropolis_single_monomer(size,p1.coordv,temperature)<<endl;
		}
	}
	//p1.outputCoordv();
	return 0;
}