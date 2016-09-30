#include "polymer.h"
#include "ranMARS.h"
#include "MonteCarlo.h"

using namespace std;
int main(){

	ploymer p1(13);
	p1.outputCoordv();
	MonteCarlo mc1;
	typedef vector<double> vd;
	vd vv(39,0);

	for(int i=0; i<vv.size()/3; i++ ){
		vv[i] = i*1.0;
		cout<<vv[i]<<endl;
	}
	vv[1]=R_Fcutoff_r;
	cout<<vv[1]<<endl;
	double Vf = sumVf(39,vv); 
	cout<<"vf:  "<<Vf<<endl;
	return 0;
}