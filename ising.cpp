#include <string>
#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <map>
#include <algorithm>
#include <sys/time.h>
//#include <stdio.h>
#include <math.h>
#include <random>

#include "Params.hpp"

#include <cmath>
#include <chrono>


//definitions
#define UP 1
#define DN -1


//Random number generator auto initialization
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> uniform_real(0.e0,1.e0) ;
std::uniform_int_distribution<> uniform_int(0,1) ;


//info regarding lattice points and spins
int L;
std::vector<std::array<int,2>> coord ;
std::map<std::array<int,2>,int> which_site;  // site index "n". ix = n%L, iy = n//L.
int nspins;
std::map<std::array<int,2>,int> spin;       // spins on the lattice

//energetics
double T, beta, J;
std::map<int,double> embetaJ ;

//observables
#define NOBS 6
#define DPQ 3   // derived physical quantities
#define ENE 0
#define ENE2 1
#define AMAG 2
#define MAG2 3
#define MAG4 4
#define MAG 5 
#define SUS 6  // susceptibility
#define CV 7   // heat capacity
#define U4 8   // binder's cumulant
std::vector<std::array<double,NOBS>> observables ;


//--------------------------------------------------------------------
void initialize_observables()
{
  observables.clear() ;
}


//--------------------------------------------------------------------
void define_Boltzmann()
{
  embetaJ.clear();

  beta = 1.e0/T;

  embetaJ[2] = exp(-2.e0*beta*J);   // required for the cluster algorithm.
  
  embetaJ[4] = exp(-4.e0*beta*J);      // we only need embeta[4] and embeta[8] because the only values that can be taken by iene is -8, -4, 0 , 4, 8 and we simply flip the spin if iene <=0, so the only boltzmann probabilities required are for iene = 4 or iene = 8.
  embetaJ[8] = exp(-8.e0*beta*J);

  //std::cout << embetaJ[4] << std::endl ;
  //std::cout << embetaJ[8] << std::endl ;
}

//--------------------------------------------------------------------	
void define_lattice_and_initialize_state()	
{	
  int is ;	
  is = 0 ;	
  for(int ix = 0; ix < L; ix++) {	
    for( int iy = 0; iy < L; iy++) {	
      coord.push_back({ix,iy}) ;	
      which_site[(std::array<int,2>){ix,iy}] = is ;	
      spin[(std::array<int,2>){ix,iy}] = UP ;	
    //spin[(std::array<int,2>){ix,iy}] = uniform_int(gen)*2 - 1;  // for a hot start with equal spin ups and downs.
      is++;	
    }	
  }	
  nspins = L*L;	
  return;	
}	
//--------------------------------------------------------------------	
void initialize_state_to_all_up()	
{	
  for(int ix = 0; ix < L; ix++) {	
    for( int iy = 0; iy < L; iy++) {	
      spin[(std::array<int,2>){ix,iy}] = UP ;	
    }	
  }	
  return;	
}  	
//--------------------------------------------------------------------	
inline int Fold(const int i, const int N) {	
  return( (i+N)%N ) ;	
}	

//--------------------------------------------------------------------
std::array<double,2> energy_density_and_magnetization()
{
  int iene=0;
  int stot=0;
  std::array<double,2> em_data;

  for(int ix = 0; ix < L; ix++) {
    for( int iy = 0; iy < L; iy++) {
      auto ixp = Fold(ix+1,L);
      auto iyp = Fold(iy+1,L);

      stot += spin[(std::array<int,2>){ix,iy}] ;
      iene += spin[(std::array<int,2>){ix,iy}]*spin[(std::array<int,2>){ixp,iy}];
      iene += spin[(std::array<int,2>){ix,iy}]*spin[(std::array<int,2>){ix,iyp}];

    }
  }

  em_data[0] = -(J*(double)iene)/((double)nspins);
  em_data[1] = ((double)stot)/((double)nspins);

  return(em_data) ;
}

//--------------------------------------------------------------------
void perform_measurements()
{
  std::array<double,2> e_m;
  e_m = energy_density_and_magnetization();

  double edens, edens2, mag, amag, mag2, mag4 ;
  edens = e_m[0];
  edens2 = edens*edens ;
  mag = e_m[1];
  amag = abs(mag);
  mag2 = mag*mag ;
  mag4 = mag2*mag2;
  observables.push_back((std::array<double,NOBS>){edens,edens2,amag,mag2,mag4,mag}) ;

  return;
}

//--------------------------------------------------------------------
// WOLFF CLUSTER ALGORITHM. (PART 1 AND 3)
//--------------------------------------------------------------------
std::vector<int> identical_nbrs(int initial_seed)   
{
  auto ix = coord[initial_seed][0] ;
  auto iy = coord[initial_seed][1] ;

  std::array<int,2> r, rpx, rmx, rpy, rmy;
  std::vector<int> id_nbrs;

  r = {ix,iy} ;
  rpx = {Fold(ix+1,L),iy} ;
  rmx = {Fold(ix-1,L),iy} ;
  rpy = {ix,Fold(iy+1,L)} ;
  rmy = {ix,Fold(iy-1,L)} ;

  if (spin[r] == spin[rpx]){
    id_nbrs.push_back(which_site[rpx]);
  }
  
  if (spin[r] == spin[rpy]){
    id_nbrs.push_back(which_site[rpy]);
  }
  
  if (spin[r] == spin[rmx]){
    id_nbrs.push_back(which_site[rmx]);
  }

  if (spin[r] == spin[rmy]){
    id_nbrs.push_back(which_site[rmy]);
  }

  return id_nbrs;
}

//--------------------------------------------------------------------
void add_to_cluster(std::vector<int> id_nbrs, std::vector<int> &C, std::vector<int> &F_new, double p)
{
  for (int j = 0; j < id_nbrs.size(); j++)
  {
    if(!(std::find(C.begin(), C.end(), id_nbrs[j]) != C.end()))   // if id_nbrs[j] doesn't belong to the vector C
    {
      if (uniform_real(gen) < p)
      {
        C.push_back(id_nbrs[j]);
        F_new.push_back(id_nbrs[j]);
      }
    }
  }  
  return;
}

//--------------------------------------------------------------------
int cluster_flip(int initial_seed)
{
  std::vector<int> C, F_old, F_new, id_nbrs;
  
  C.push_back(initial_seed);
  F_old.push_back(initial_seed);

  double p = 1 - embetaJ[2];

  while (F_old.size()!=0)
  {
    F_new.clear();
    for (int i = 0; i < F_old.size(); i++)
    {
      id_nbrs = identical_nbrs(F_old[i]);
      add_to_cluster(id_nbrs, C, F_new, p);       // updates C and F_new arrays after selecting new spins into the cluster.
    }
    F_old = F_new;
  }

  for (int i = 0; i < C.size(); i ++ )
  {
    int ix = (int)(C[i]/L);
    int iy = (int)(C[i]%L);
    spin[(std::array<int,2>){ix,iy}]*=-1;   // spin flip of the cluster
  }

  return C.size();
}

//--------------------------------------------------------------------
double run_wolff_cluster(int const neqsweeps, int const nsamsweeps, int const nsampstp)
{
  std::vector<int> C_sizes;     // vector of cluster sizes
  int C_len;             

  // Equilibration process. 
  for(int isw = 0 ; isw < neqsweeps ; isw++)
  {
    auto initial_seed = (int)(uniform_real(gen)*(double)nspins);
    C_len = cluster_flip(initial_seed);
//  C_sizes.push_back(C_len);
  }
 
 // Thermalized.
  for(int isw = 0 ; isw < nsamsweeps ; isw++)
  {
    auto initial_seed = (int)(uniform_real(gen)*(double)nspins);
    C_len = cluster_flip(initial_seed);
    C_sizes.push_back(C_len);                   // collecting cluster sizes for equilibrated sweeps.

    if( isw%nsampstp == 0 ) perform_measurements() ;
  }

  double C_avg = 0.0e0;
  for (int i = 0; i < C_sizes.size(); i++)
  {
    C_avg += C_sizes[i]; 
  }
  C_avg /= C_sizes.size();

  return C_avg;             // average cluster size in the MC run (required for appropriate normalization of autocorr_time).
}


//--------------------------------------------------------------------
// METROPOLIS-HASTINGS ALGORITHM. (PART 1 AND 3)
//--------------------------------------------------------------------
void flip_spin(int const ix, int const iy)
{
  int iene = 0;

  std::array<int,2> r, rpx, rmx, rpy, rmy;

  r = {ix,iy} ;
  rpx = {Fold(ix+1,L),iy} ;
  rmx = {Fold(ix-1,L),iy} ;
  rpy = {ix,Fold(iy+1,L)} ;
  rmy = {ix,Fold(iy-1,L)} ;

  iene = 2*spin[r]*(spin[rpx]+spin[rmx]+spin[rpy]+spin[rmy]);


  if ( (iene <= 0) || ( uniform_real(gen) < embetaJ[iene]) ) spin[r] *= -1;
  return;
}

//--------------------------------------------------------------------
	void monte_carlo_sweep_sequential()	
{	
  for(int ix = 0; ix < L; ix++) {	
    for( int iy = 0; iy < L; iy++) {	
      flip_spin(ix,iy);	
    }	
  }	
  return ;	
}

//--------------------------------------------------------------------
void monte_carlo_sweep_random()
{

  for(int is = 0; is < nspins ; is++) {
    auto spin_to_flip = (int)(uniform_real(gen)*(double)nspins);
    auto ix = coord[spin_to_flip][0] ;
    auto iy = coord[spin_to_flip][1] ;
    flip_spin(ix,iy);
  }

  return ;
}

//--------------------------------------------------------------------
void run_metropolis(int const neqsweeps, int const nsamsweeps, int const nsampstp, bool const sequential)
{
  //Equilibration
  //std::cout << "Running " << neqsweeps << " equilibrium sweeps..." ;
  for(int isw = 0 ; isw < neqsweeps ; isw++){
    if(sequential) {
      monte_carlo_sweep_sequential() ;
    }
    else{
      monte_carlo_sweep_random() ;
    }
      
  }
  //std::cout << "...Done. " << std::endl ;

  //std::cout << "Running " << nsamsweeps << " sampling sweeps..." ;
  for(int isw = 0 ; isw < nsamsweeps ; isw++){
    if(sequential) {
      monte_carlo_sweep_sequential() ;
    }
    else{
      monte_carlo_sweep_random() ;
    }
    if( isw%nsampstp == 0 ) perform_measurements() ;
  }
  //std::cout << "...Done. " << std::endl ;
  
}

//--------------------------------------------------------------------
// AUTOCORRELATION MEASUREMENTS. (PART 2)
//--------------------------------------------------------------------
double obtain_correlation_value(int iobs, int t)
{
  double mean = 0.e0, mean2 = 0.e0;
  double val;
  double corsum = 0.e0;
  for (int k = 0; k < observables.size()-t; k++)    // here.
  {
      val = observables[k][iobs];
      mean += val;
      mean2 += val*val;
  }
  mean /= (double)(observables.size()-t);    // here.
  mean2 /= (double)(observables.size()-t);   // here.

  for (int j = 0; j < observables.size()-t; j++)
  {
      corsum += (observables[j][iobs] - mean)*(observables[j+t][iobs] - mean);
  }
  double C = (corsum/((observables.size()-t)*(mean2 - mean*mean)));  // correlation function value for time lag t.
  return C;
}
//--------------------------------------------------------------------

double obtain_correlation_time(int const tmax, std::ofstream &fp, std::ofstream &fo)
{
  // saving observations in a data file
  for (int it = 0; it < observables.size(); it ++)            // fo is the new ofstream variable that I'm using to output a file of observable measurements.
  {
    fo << it << "," << T << "," ;
    for (int iobs = 0; iobs < NOBS ; iobs ++)
    {
      fo << observables[it][iobs] << ",";
    }
    fo << std::endl;
  }

  auto ndat = observables.size() - tmax;

  if (ndat <= 0 ) {
    std::cout << observables.size() << " < tmax " << tmax << " cannot compute correlations" ;
    return 0;
  }
  
  std::vector<std::vector<double>> cor;

  for (int iobs = 0; iobs < NOBS - 1; iobs++)
  {
      std::vector<double> ct;
      for (int t = 0; t <= tmax; t++)
      {
          double C = obtain_correlation_value(iobs, t);
          ct.push_back(C);
      }
      cor.push_back(ct);
  }

  fp << std::scientific ; 
  for(int t=0; t <= tmax; t++) 
  {
    fp << t << "," << T << "," ;
    for(int iobs = 0; iobs < NOBS - 1; iobs++) 
    {
      fp << cor[iobs][t] << "," ;
    }
    fp << std::endl ;
  }

  fp << std::endl ;

  // Calculate the Integrated Auto-Correlation time using "automatic windowing" procedure. (Ref. A. Sokal Monte Carlo Methods in Statistical Mechanics: Foundations and New Algorithms  https://link.springer.com/chapter/10.1007/978-1-4899-0319-8_6)

  std::array<double,NOBS-1> tau_array ;  // store autocorrelation times for each observables at a given (T, L) (Except magnetization)
  
  for(int iobs = 0; iobs < NOBS-1; iobs++) {
    double tau_sum = 0.50e0 ; 
    int N = 20;                    // start with some random cut-off and keep changing it accordingly if it's smaller than 6*tau.
    bool condition = true;
    int it = 1;

    while (condition)
    {
      for(it; it < N; it++) {  
        tau_sum = tau_sum + (cor[iobs][it]);
      }
      condition = N <= (6*tau_sum);
      // update N.
      N = N*1.3;
      if (N > tmax) 
      {
        std::cout <<std::endl<<std::endl << "N > tmax ("<< N << ">" << tmax <<"). Choose a larger value of tmax" << std::endl ;
        //exit(1);              // edit this out later.
        return tau_sum;
      }
    }
    N = N/1.2; //undoing the last update.
    tau_array[iobs] = tau_sum;

    std::cout << "tau_sum = " << tau_sum << "  N = " << N << std::endl; // comment it out later, just for troubleshoot.
  }

  double tau_int = *std::max_element(tau_array.begin(), tau_array.end()); // maxima of correlation time out of all observables

  return tau_int;  // return the integrated correlation time
}

//--------------------------------------------------------------------
// DATA ANALYSIS. (PART 4)
//--------------------------------------------------------------------
double mean_data(std::vector<double> data)
{
  double average = 0;
  for (int i = 0; i < data.size(); i++)
  {
    average += data[i];
  }
  average/=data.size();

  return average;
}

//--------------------------------------------------------------------
double jackknife_error(std::vector<double> data, double mean)
{
  double summ = 0;
  for (int i = 0; i < data.size(); i++)
  {
    summ += (data[i] - mean)*(data[i] - mean); 
  }

  double fact = (double)(data.size()-1)/(data.size());
  double jack_error = sqrt(fact*summ);
  
  return jack_error;
}

//--------------------------------------------------------------------
void analyze_samples(int const nblen, std::ofstream &fp, double cortime, std::ofstream &fu)
{
  // saving uncorrelated measurements in a data file
  for (int it = 0; it < observables.size(); it ++)     
  {
    fu << it << "," << T << "," ;
    for (int iobs = 0; iobs < NOBS ; iobs ++)
    {
      fu << observables[it][iobs] << ",";
    }
    fu << std::endl;
  }
  
  auto ndat = observables.size() ; 
  
  std::vector<std::vector<double>> bin_avg;   
  std::vector <double> obs_bin_avg;
    
  for (int iobs = 0; iobs < NOBS; iobs++)
  {
    double avg = 0;
    obs_bin_avg.clear();
    double count = 0;
    for (int i = 0; i < ndat; i++)
    {
      avg += observables[i][iobs];
      count +=1;
            
      if ((i+1)%nblen == 0 || i == ndat - 1)
      {
        avg /= count;
        obs_bin_avg.push_back(avg);
        avg = 0;
        count = 0;
      }
    }
  bin_avg.push_back(obs_bin_avg);
  }

  double nbins = bin_avg[0].size();
  // bin_avg[iobs][j] is the mean of elements in bin "j" of observable "iobs" measurements.

  // Now we also need to insert the binned averages of derived physical quantities (chi, C_v, U_4)into this vector.

  std::vector<double> susc;
  std::vector<double> heatcap;
  std::vector<double> bindercum;

  for (int j = 0; j < nbins; j++)
  {
    susc.push_back(beta*(double)nspins*(bin_avg[MAG2][j] - bin_avg[AMAG][j]*bin_avg[AMAG][j]));
    heatcap.push_back(beta*beta*(double)nspins*(bin_avg[ENE2][j] - bin_avg[ENE][j]*bin_avg[ENE][j]));
    bindercum.push_back(1.e0-(bin_avg[MAG4][j]/(3.e0*bin_avg[MAG2][j]*bin_avg[MAG2][j])));
  }

  bin_avg.push_back(susc);
  bin_avg.push_back(heatcap);
  bin_avg.push_back(bindercum);

  // At this point, bin_avg[iobs][j] is the mean of elements in bin "j" of the physical quantity "iobs" measurements. The label for indexing physical quantities i.e. "iobs" goes from 0 to 8, with 6, 7, 8 now representing susceptibility, heat capacity, and binder's cumulant.

  std::array<double, NOBS+DPQ> obs_dpq_means;
  for (int iobs=0; iobs < NOBS+DPQ; iobs++)
  {
    obs_dpq_means[iobs] = mean_data(bin_avg[iobs]);
  }
  // Now we also have the grand means of all the observables + derived quantities.

  
  std::vector<std::vector<double>> jackknife_bin_avg;
  std::vector<double> temp;
  for (int iobs = 0; iobs < NOBS+DPQ; iobs++)
  {
    temp.clear();
    for (int j = 0; j < nbins; j++)
    {
      temp.push_back( (ndat*obs_dpq_means[iobs] - nblen*bin_avg[iobs][j])/(ndat - nblen) );
    }
    jackknife_bin_avg.push_back(temp);
  }

  std::array<double, NOBS+DPQ> obs_dpq_jackknife_errors;
  for (int iobs = 0; iobs < NOBS+DPQ; iobs++)
  {
    obs_dpq_jackknife_errors[iobs] = jackknife_error(jackknife_bin_avg[iobs], obs_dpq_means[iobs]);
  }

  // Printing and writing the means and the errors obtained using JackKnife.
  std::cout << T << " " << L << " " << cortime << " "; 
  std::cout << observables.size() << " " << nblen << " " << nbins << " ";
  std::cout << obs_dpq_means[ENE] << " " << obs_dpq_jackknife_errors[ENE] << " " ;
  std::cout << obs_dpq_means[ENE2] << " " << obs_dpq_jackknife_errors[ENE2] << " " ;
  std::cout << obs_dpq_means[AMAG] << " " << obs_dpq_jackknife_errors[AMAG] << " " ;
  std::cout << obs_dpq_means[MAG2] << " " << obs_dpq_jackknife_errors[MAG2] << " " ;
  std::cout << obs_dpq_means[MAG4] << " " << obs_dpq_jackknife_errors[MAG4] << " " ;
  std::cout << obs_dpq_means[MAG] << " " << obs_dpq_jackknife_errors[MAG] << " " ;
  std::cout << obs_dpq_means[SUS] << " " << obs_dpq_jackknife_errors[SUS] << " " ;
  std::cout << obs_dpq_means[CV] << " " << obs_dpq_jackknife_errors[CV] << " " ;
  std::cout << obs_dpq_means[U4] << " " << obs_dpq_jackknife_errors[U4] << " " ;
  std::cout << std::endl;


  fp << std::scientific ;
  fp << T << "," << L << "," << cortime << ","; 
  fp << observables.size() << "," << nblen << "," << nbins << ",";
  fp << obs_dpq_means[ENE] << "," << obs_dpq_jackknife_errors[ENE] << "," ;
  fp << obs_dpq_means[ENE2] << "," << obs_dpq_jackknife_errors[ENE2] << "," ;
  fp << obs_dpq_means[AMAG] << "," << obs_dpq_jackknife_errors[AMAG] << "," ;
  fp << obs_dpq_means[MAG2] << "," << obs_dpq_jackknife_errors[MAG2] << "," ;
  fp << obs_dpq_means[MAG4] << "," << obs_dpq_jackknife_errors[MAG4] << "," ;
  fp << obs_dpq_means[MAG] << "," << obs_dpq_jackknife_errors[MAG] << "," ;
  fp << obs_dpq_means[SUS] << "," << obs_dpq_jackknife_errors[SUS] << "," ;
  fp << obs_dpq_means[CV] << "," << obs_dpq_jackknife_errors[CV] << "," ;
  fp << obs_dpq_means[U4] << "," << obs_dpq_jackknife_errors[U4] << "," ;
  fp << std::endl;

  return;
}

//--------------------------------------------------------------------
// MAIN FUNCTION.
//--------------------------------------------------------------------

int main(int argc, char** argv)
{
  //PramsReader reader;
  Params parameters;
  
  std::string casename;

  casename = "tis" ;
  
  //lattice
  L=64;

  //energetics
  J=1.e0;
  // T=2.2e0;
  // beta=1.e0/T;

  //montecarlo
  int neqsweeps, nsamsweeps, nsampstp, nblen;
  neqsweeps = 512;
  nsamsweeps = 32768;
  nsampstp = 32 ;
  nblen = 16;

  bool sequential = false;
  bool clust = true;
  
  double Tmin=2.1e0, Tmax=2.40e0, dT=0.02e0;

  bool initialize_all_up = false, print_cor = false;
  int tmax=5000;
  int N_out = 500;  // number of measurements desired while using stepsize = 2*tau.
  // NOT READING N_out from the in.ising.json file for some reason.
  
  std::cin >> parameters ;

  std::cout << parameters.toStyledString() << std::endl ;
  
  //Initialize random seed
  //_SeedRandom() ;
  
  for(auto card = parameters.begin(); card != parameters.end(); ++card) {
    auto cardname = card.key();
    if(cardname == "casename") {
      std::cout << "found casename" << std::endl;
      casename = parameters["casename"].asString() ;
      std::cout << "Running case : " << casename << std::endl ; 
    }
    else if(cardname == "L") {
      std::cout << "found L" << std::endl ;
      L = parameters["L"].asInt() ;
      std::cout << "L = " << L << std::endl ;
    }
    else if( cardname == "Tmin" ) {
     std::cout << "found Tmin" << std::endl ;
     Tmin = parameters["Tmin"].asDouble() ;
     std::cout << "Tmin = " << Tmin << std::endl ;
    }
    else if( cardname == "Tmax" ) {
      std::cout << "found Tmax" << std::endl ;
      Tmax = parameters["Tmax"].asDouble() ;
      std::cout << "Tmax = " << Tmax << std::endl ;
    }
    else if( cardname == "dT" ) {
      std::cout << "found dT" << std::endl ;
      dT = parameters["dT"].asDouble() ;
      std::cout << "dT = " << dT << std::endl ;
    }
    else if( cardname == "neqsweeps" ) {
     std::cout << "found neqsweeps" << std::endl ;
     neqsweeps = parameters["neqsweeps"].asInt() ;
     std::cout << "neqsweeps = " << neqsweeps << std::endl ;
    }
    else if( cardname == "nsamsweeps" ) {
     std::cout << "found nsamsweeps" << std::endl ;
     nsamsweeps = parameters["nsamsweeps"].asInt() ;
     std::cout << "nsamsweeps = " << nsamsweeps << std::endl ;
    }
    else if( cardname == "nsampstp" ) {
     std::cout << "found nsampstp" << std::endl ;
     nsampstp = parameters["nsampstp"].asInt() ;
     std::cout << "nsampstp = " << nsampstp << std::endl ;
    }
    else if( cardname == "nblen" ) {
      std::cout << "found nblen" << std::endl ;
      nblen = parameters["nblen"].asInt() ;
      std::cout << "nblen = " << nblen << std::endl ;
    }
    else if( cardname == "initialize_all_up" ) {
      std::cout << "found initialize_all_up" << std::endl ;
      initialize_all_up = parameters["initialize_all_up"].asBool() ;
      std::cout << "initialize_all_up = " << initialize_all_up << std::endl ;
    }
    else if( cardname == "print_cor" ) {
      std::cout << "found print_cor" << std::endl ;
      print_cor = parameters["print_cor"].asBool() ;
      std::cout << "print_cor = " << print_cor << std::endl ;
    }
    else if( cardname == "sequential" ) {
      std::cout << "found sequential" << std::endl ;
      sequential = parameters["sequential"].asBool() ;
      std::cout << "sequential = " << sequential << std::endl ;
    }
    else if( cardname == "cluster" ) {
      std::cout << "found cluster" << std::endl ;
      clust = parameters["cluster"].asBool() ;
      std::cout << "cluster = " << clust << std::endl ;
    }
    else if( cardname == "tmax" ) {
      std::cout << "found tmax" << std::endl ;
      tmax = parameters["tmax"].asInt() ;
      std::cout << "tmax = " << tmax << std::endl ;
    }
    else if( cardname == "N_out" ) {
      std::cout << "found N_out" << std::endl ;
      N_out   = parameters["N_out"].asInt() ;
      std::cout << "N_out = " << N_out << std::endl ;
    }
    else {
      std::cout << "Unknown parameter " << cardname << std::endl;
      exit(-1) ;	
    }
   
  }
  
  //Define Lattice
  define_lattice_and_initialize_state() ;

  //createoputfile
  auto filename = "output/"+casename+"_out.plt" ;
  std::ofstream foutp(filename);

  auto measurement_file ="output/"+casename+"_measurements.plt" ;
  std::ofstream fobs(measurement_file);

  auto uncor_data = "output/" + casename + "_uncorrdata.plt";
  std::ofstream funcorr(uncor_data);

  //if(print_cor) {
    auto corfile = "output/"+casename+"_cor.plt" ;
    std::ofstream fcor(corfile) ;
    //}

  T = Tmin;
  
  while (T <= Tmax+dT/2.e0) { 
    double C_avg = 0.0e0;

    std::cout << "--------------------------------------------------------------" << std::endl;
    std::cout << "L = " << L << ", T = " << T << std::endl;
    std::cout << "PART 1" << std::endl;
    // preliminary monte carlo (part 1).

    auto start1 = std::chrono::steady_clock::now();
    define_Boltzmann();
    initialize_observables();
    if(initialize_all_up) initialize_state_to_all_up() ;    

    if (clust)    {
      C_avg = run_wolff_cluster(neqsweeps, nsamsweeps, 1);
    }
    else    {
      run_metropolis(neqsweeps, nsamsweeps, 1, sequential) ;
    }
    auto end1 = std::chrono::steady_clock::now();

    std::cout<<std::chrono::duration_cast<std::chrono::seconds>(end1-start1).count()<< " seconds for Monte Carlo part 1 to calculate tau.  (sweeps = " << nsamsweeps << ")" << std::endl<< std::endl;


    //find and print autocorrleations (part 2).
    std::cout << "PART 2" << std::endl;
    auto start2 = std::chrono::steady_clock::now();
    double tau;
    tau = obtain_correlation_time(tmax, fcor, fobs) ;
    if (clust)    {
      std::cout << "<c> (avg cluster size) = " << C_avg << std::endl;
      std::cout << "τ (scaled by <c>/N_spins) = " << (double)(tau*(C_avg/nspins)) << std::endl ;
     // rescaling the time unit for wolff cluster for relevant comparisions of times.
    }
    else    {
      std::cout << "τ = " << tau << std::endl ;
    }
    auto end2 = std::chrono::steady_clock::now();

    std::cout<<std::chrono::duration_cast<std::chrono::seconds>(end2-start2).count()<< " seconds for calculating autocorrelation time τ. " << std::endl<< std::endl;

    // running monte carlo with known autocorrelation times (part 3).
    std::cout << "PART 3" << std::endl;
    define_Boltzmann();
    initialize_observables();
    if(initialize_all_up) initialize_state_to_all_up() ;    

    int nsampstp_new = 2*tau;
    int nsamsweeps_new = nsampstp_new*N_out;
    auto start3 = std::chrono::steady_clock::now(); 
    if (clust)    {
      C_avg = run_wolff_cluster(neqsweeps, nsamsweeps_new, nsampstp_new);
    }
    else    {
      run_metropolis(neqsweeps, nsamsweeps_new, nsampstp_new, sequential);
    }       
    auto end3 = std::chrono::steady_clock::now();

    std::cout<<std::chrono::duration_cast<std::chrono::seconds>(end3-start3).count()<< " seconds for the Monte Carlo part 2 to calculate expvals and errors. (sweeps = " << nsamsweeps_new << " = " << nsampstp_new << "*" << N_out << ")" <<std::endl << std::endl;

    // data analysis (part 4).
    std::cout << "PART 4" << std::endl;
    auto start4 = std::chrono::steady_clock::now();    
    if (clust)    {
      analyze_samples(nblen, foutp, (double)(tau*(C_avg/nspins)), funcorr) ;
    }
    else    {
      analyze_samples(nblen, foutp, tau, funcorr);
    }
    auto end4 = std::chrono::steady_clock::now();

    std::cout<<std::chrono::duration_cast<std::chrono::seconds>(end4-start4).count()<< " seconds for data analysis. " << std::endl << std::endl<< std::endl;

    T+=dT;
  }
  fobs.close() ;
  foutp.close() ;
  fcor.close() ;   
  return(1) ;
}
  