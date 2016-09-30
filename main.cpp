#include "polymer.h"
#include "ranMARS.h"
#include "MonteCarlo.h"

using namespace std;
int main(){

	ploymer p1(13);
	p1.outputCoordv();
	MonteCarlo mc1;

	typedef vector<double> vd;
	vd vv(40,0);

	for(int i=0; i<vv.size()/3; i++ ){
		vv[i] = i*1.0;
		cout<<vv[i]<<endl;
	}
	vv[1]=1.2;
	cout<<vv[1]<<endl;
	int size = floor((vv.size()+1)/3);
	cout<<"size is "<<size<<endl;

	double Vf = sumVf(size,vv); 
	cout<<"vf:  "<<setw(15)<< setprecision(10)<<Vf<<endl;
	sumVlj(6,vv);
	double V_lj=0.0, V_fene = 1.0;
	cout <<"V_lj: "<<V_lj <<endl;
	return 0;
}