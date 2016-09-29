/* the class is to create a polymers class the for Monte Carlo update
1) geo-properties
-- n : number of mers
-- coordnations[3*n+1]: 3n for (x,y,z)*n, 1 for possible index
2) thermo properties
-- T  : temperature
-- E  : total E
-- Cv : the specific heat
-- Rg :
*/


# include <vector>
# include "ranMARS.h" 
# include <iostream>
# include <iomanip> // setw()

using namespace std;
class ploymer{
public: 
	int n =0;
	vector<double> coordv; // coordinate vector
	double T=0.0;
	double E=0.0;
// constructor1, give the n value and set the coordations as (0,0,0), (1,0,0) ,~~~(n,0,0);
	ploymer(int number);
// constructor2, give the n value and then give the coordinates in a file
//	ploymer(double n, filepointer *file);


// metropolis updata your system by change 1 mer's position 
	double metropolis_1(double Eold, double T);
// metropolis updata your system by change n mer's position 
	double metropolis_n(double Eold, double T, int n);
// ouput coordinate
	void outputCoordv();

};

ploymer::ploymer(int number){
	if (number<=0) { 
		cout<<"Eeror ploymer's number <= 0\n"; 
		exit(1);
	}
	n = number;
	// coordinator
	coordv.assign(3*n+1,0.0);
	for (int i=0; i<n; i++){
		coordv[i]=i;						// x axis only, coordv-> (x1,x2~~),(y1,y2,~~),(z1,z2~~)
		//cout <<i<<"th "<<coordv[i]<<endl; 
	}
}


void ploymer::outputCoordv(){
	for (int i=0; i<n; i++){
		cout <<setw(3)<<i+1
		<<"th mer ("
		<<setw(3)<<coordv[i]<<" , "
		<<setw(3)<<coordv[i+n]<<" , "
		<<setw(3)<<coordv[i+n+n]<<")"<<endl;
	}
}
		 
