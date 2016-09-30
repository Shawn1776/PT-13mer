/**	Author: Shawn Zhang @uga physics 
  *	Date  : Sep 2016
  */
#include "ranMARS.h"
#include <vector>
#include <math.h>
#include <limits.h>

// for Fene; maximum deviation is r-ro = R
#define K (98.0/5.0) 					// spring constant, in(energy/distance^2)
#define R (3.0/7.0)  					//  be careful of the BRACKE!!!  // (distance)
#define r0 1.0     						// equilibrium bond length (distance)
#define R_Fcutoff_l 0.571428571428572   // the left boundary of the FENE potential
#define R_Fcutoff_r 1.428571428571428   // FENE right boundary  
#define DBL_MAX 1.7976931348623158e+308 // max value in double type data, return this value is r is outoff above range 

// constant for LJ potential
#define EPSILON 1.0                     // epsilon (energy)
#define SIGMA (1.0/(pow(2,1.0/6.0)))    // sigma (distance)    (1.0/1.122462048309373 = 0.8908987181403393)
#define epsilon_LJ 1.0
#define V_LJ_cut -0.016316891136        // energy cutoff for LJ potential, accurate to the 12th digit,,,
#define R_LjCutoff_l 0.1                // LJ left boundary, it will return a huge U_Lj energy number, but it is toke care by the energy boundray dE/T < 600
#define R_LjCutoff_r 2.22724679535085   // LJ right boundary, a cutting off distance which we choose to make the Lennard-Jones potential go to 0 when the distance between two beads is larger than r_c


using namespace std;

double calcdist(double a, double b, double c, double x, double y, double z);
// /** calc_dist_v
//   * 
//   */
// std::vector<double> calc_dist_v();
double  sumVf(int size, std::vector<double> &v);
double sumVlj(int size, std::vector<double> &v);// input the a[][] Or b[][]


class MonteCarlo{
 	public:
 	// Metroplis update a single monomer of a polymer
	double Metroplis_single_monomer(vector<double> &v, double T, int seed);


 };


// Metroplis_single_monomer
 // precondition: 1) structure-> calculate the totalE = Vlj+Vfene
 //				  2) temperature
//				  3) random number
double MonteCarlo::Metroplis_single_monomer(vector<double> &v, double T, double random){
	double V_lj=0.0, V_fene = 0.0, Eold =0.0, Enew =0.0;

	cout <<"V_lj: "<<V_lj <<endl;

	return 0.0;
}
 
double calcdist( double a, double b, double c, double x, double y, double z){	
   return sqrt(pow(a - x, 2) + pow(b - y, 2) + pow(c - z, 2));
}

// calclate the sum of Vfene on given the coord vector
double sumVf(int size, std::vector<double> &v){
	double Vfene = 0.0;
	double r = 0.0; 
	double coefficient = -0.5*K*R*R;		// coefficient = 1/2.0 *K*R^2
	for (int i=0; i<(size-1); i++){
		r = calcdist(v[i],v[i+size],v[i+size+size],v[i+1],v[i+1+size],v[i+1+size+size]);
		if(r <= R_Fcutoff_l || r >= R_Fcutoff_r){
			return DBL_MAX;
		}else{
			Vfene += coefficient*log(1.0-(r-r0)*(r-r0)/R/R);
			// cout<<i<<" th distance: "<<setw(15)<< setprecision(10)<<Vfene<<" r is "<<r<<endl; // debug
		}//if elif
	}//for 
		cout <<"sigma: "<<SIGMA<<endl;
		return Vfene;	
}


// calclate the sum of V_Lj on given the coord vector
double sumVlj(int size, std::vector<double> &v){
	double VLj =0.0,r=0.0;
	int cnt =0;  // i,j when j is not i+1 or i-1
	for (int i = 0; i<size-2; i++){
		for(int j = i+2; j<size; j++){
			cnt++;
			r = calcdist(v[i],v[i+size],v[i+size+size],v[j],v[j+size],v[j+size+size]);
			if(r <= R_LjCutoff_l){ // if r is smaller than left boundary, then, return max energy
				return DBL_MAX;
			}else if(r <= R_LjCutoff_r){ // if r exceeds the right boundary then, return 0 energy
				VLj += 4.0*EPSILON*(pow((SIGMA/r),12.0)-pow((SIGMA/r),6.0))-V_LJ_cut;
			} // if -else if
		}//j
	}//i
	return VLj;
}