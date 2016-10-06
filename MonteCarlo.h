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
// calculate V_fene
double  sumVf(int size, std::vector<double> &v); // pass size to optimize the code, we donot have to since we can calculate from the vector
// calculate V_Lj
double sumVlj(int size, std::vector<double> &v); // pass size to optimize the code
// check if V_fene or V_Lj are still in boundary
bool is_DBL_MAX(double v_fene, double V_Lj);


class MonteCarlo{
 	public:
 	// Metroplis update a single monomer of a polymer
	double Metropolis_single_monomer(int size, vector<double> &v, double T);  // pass size to optimize the code
 };

// method: Metroplis_single_monomer
//  precondition: 1) structure-> calculate the totalE = Vlj+Vfene
//				  2) temperature
//				  3) random number ( can generated from main function then pass to this function)
//  postcodition: if update accepted, return new E; else return old E 
double MonteCarlo::Metropolis_single_monomer(int size, vector<double> &v, double T){ 	
	double V_lj=0.0, V_fene = 0.0, Eold =0.0, Enew =0.0;
	Eold = sumVf(size, v)+sumVlj(size,v);
	// random select 1 atom from 13 atoms
	int index = floor(RAN01()*13);                // initalize at the main function
	// update x,y,z of index atom, first save the old coords
	double coord_old[3] = {}; // inital the array to zeros
	for (int i=0; i<3; i++){
		// first save the old coords
		coord_old[i] = v[index+13*i];
		// generate the new ones within a certain range (=/-0.3 in this case)
		v[index+13*i] += 0.3*(RAN01()-0.5);
	}
	// calculation new E
	double Vfene_new=sumVf(size, v), Vlj_new = sumVlj(size,v); // could be optimized on the V_fene()
	if (!is_DBL_MAX(Vfene_new,Vlj_new)){
		// cout <<"do work"<<endl;
		// if V_fene and V_Lj is not out off boundary then do~~~
		Enew = sumVf(size, v)+sumVlj(size,v);
	}else{
		//if V_fene and V_Lj is not out off boundary then re-roll system
// cout<<"debug V_fene or V_Lj out off boundary"
		for (int i=0; i<3; i++){
				v[index+13*i]=coord_old[i];
			}
			// return Eold
			return Eold;
	}
	// metropolis update according to the Metropolis criterion
	double deltaE = Enew-Eold;
	if (deltaE<=0.0){
		return Enew;
	}else{
		double p_acc = exp(-1.0*deltaE/T);    // accepted probobility P_acc
//		cout <<"p_acct "<<p_acc<<endl;
		double r = RAN01();                   // random number
		//r=0.1;
		if( p_acc>r ){ // accepted 
			return Enew;
		}else{        // denied 
			//system rolls back
			for (int i=0; i<3; i++){
				v[index+13*i]=coord_old[i];
			}
			// return Eold
			return Eold;
		} // iner if-else
	}// outer if-else
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
//		cout <<"sigma: "<<SIGMA<<endl;
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


// check if the Vlj or Vfene is DBL_MAX, if is return ture, else return false
bool is_DBL_MAX(double v_fene, double V_Lj){
	if (DBL_MAX == v_fene) 
		return true;
	if (DBL_MAX == V_Lj)
		return true;
	if (DBL_MAX <= v_fene+V_Lj) // check if the sum is outoff boundary, Oct062016, ShawnZ
		return true;
	return false;
}