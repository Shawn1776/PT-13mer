/**	Author: Shawn Zhang @uga physics 
  *	Date  : Sep 2016
  */
#include "ranMARS.h"
#include <vector>
#include <math.h>
#include <limits.h>


// constant for MC LJ and FENE ///////////////////////////////////////////
#define EPSILON 1 // epsilon is the constant for Lernord jones
// Sigma is constant for Lernord Jones// which is         1.0/1.12246205 = 0.832490
#define SIGMA 1.0/(pow(2,1.0/6.0)) 

#define K 98.0/5.0
// for Fene; maximum deviation is r-ro = R
#define R 3.0/7.0

#define r0 1.0
// max value in double type data
#define DBL_MAX 1.7976931348623158e+308

#define epsilon_LJ 1.0
#define V_LJ_cut -0.016316891136

//LJ left boundary // it is a huge number, but it is toke care by the energy boundray dE/T < 600
#define R_LjCutoff_l 1.0                   
// LJ right boundary
#define R_LjCutoff_r 2.22724679535085 //A cutting off distance which we choose to make the Lennard-Jones potential go to 0 when the distance between two beads is larger than r_c
#define r_LJ_cut 0.1 //artificial cut off of LJ potential 

// the left boundary of the FENE potential
#define R_Fcutoff_l 0.571428571428572   
// FENE right boundary  
#define R_Fcutoff_r 1.428571428571428

// constant for MC LJ and FENE ///////////////////////////////////////////

using namespace std;
double calcdist(double a, double b, double c, double x, double y, double z);

// /** calc_dist_v
//   * 
//   */
// std::vector<double> calc_dist_v();

double sumVf( int n, std::vector<double> &v);
double sumVlj(int n, std::vector<double> &v);// input the a[][] Or b[][]




 class MonteCarlo{


 public:
 	// Metroplis update a single monomer of a polymer
	double Metroplis_single_monomer(vector<double> &v, double T, int seed);

 };


 
double MonteCarlo::Metroplis_single_monomer(vector<double> &v, double T, int seed){


}




// 
double calcdist( double a, double b, double c, double x, double y, double z){	
   return sqrt(pow(a - x, 2) + pow(b - y, 2) + pow(c - z, 2));
}

// calclate the sum of Vfene on given the coord vector
double sumVf( int n, std::vector<double> &v){
	double Vfene = 0.0;
	int size = (v.size()+1)/3;
	cout<<"size is "<<size<<endl;
	double r = 0.0; 
	double coefficient = (-1.0/2.0)*K*pow(R,2);
	for (int i=0; i<(size-1); i++){
		r = calcdist(v[i],v[i+size],v[i+size+size],v[i+1],v[i+1+size],v[i+1+size+size]);
		cout <<"r "<<r<<endl;
		// r is in between (R_Fcutoff_l and R_Fcutoff_r)
		if(r<=R_Fcutoff_l || r>=R_Fcutoff_r){
			return DBL_MAX;
		}else{
			Vfene = Vfene+log(1.0-pow((r-r0)/(R),2));
		}//if elif
	}//for 
		return Vfene;	
}


// calclate the sum of V_Lj on given the coord vector
double sumVlj(int n, std::vector<double> &v){
	return 0.0;
}