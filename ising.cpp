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


//definitions
#define UP 1
#define DN -1


//Random number generator auto initialization
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> uniform(0.e0,1.e0) ;



//info regarding lattice points and spins
int L;
std::vector<std::array<int,2>> coord ;
std::map<std::array<int,2>,int> which_site;
int nspins;
std::map<std::array<int,2>,int> spin;

//energetics
double T, beta, J;
std::map<int,double> embetaJ ;

//observables
#define NOBS 6
#define ENE 0
#define ENE2 1
#define MAG 2
#define AMAG 3
#define MAG2 4
#define MAG4 5
std::vector<std::array<double,NOBS>> observables ;


//--------------------------------------------------------------------
void initialize_observables()
{
  observables.clear() ;
}


//--------------------------------------------------------------------
void define_Boltzmann( )
{
  embetaJ.clear() ;

  beta = 1.e0/T;
  
  embetaJ[4] = exp(-4.e0*beta*J);
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
double energy_density()
{
  int iene=0;

  for(int ix = 0; ix < L; ix++) {
    for( int iy = 0; iy < L; iy++) {
      auto ixp = Fold(ix+1,L);
      auto iyp = Fold(iy+1,L);
      //xbond
      iene += spin[(std::array<int,2>){ix,iy}]*spin[(std::array<int,2>){ixp,iy}];
      //ybond
      iene += spin[(std::array<int,2>){ix,iy}]*spin[(std::array<int,2>){ix,iyp}];
    }
  }

  return(-(J*(double)iene)/((double)nspins)) ;
}


//--------------------------------------------------------------------
double magnetization()
{
  int stot=0;

  for(int ix = 0; ix < L; ix++) {
    for( int iy = 0; iy < L; iy++) {
      auto ixp = Fold(ix+1,L);
      auto iyp = Fold(iy+1,L);
      //xbond
      stot += spin[(std::array<int,2>){ix,iy}] ;
    }
  }

  return(((double)stot)/((double)nspins)) ;
}

//--------------------------------------------------------------------
void perform_measurements()
{
  double edens, edens2, mag, amag, mag2, mag4 ;
  edens = energy_density();
  edens2 = edens*edens ;
  mag = magnetization();
  amag = abs(mag);
  mag2 = mag*mag ;
  mag4 = mag2*mag2;
  observables.push_back((std::array<double,NOBS>){edens,edens2,mag,amag,mag2,mag4}) ;

  return;
}

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

  
  if ( (iene <= 0) || ( uniform(gen) < embetaJ[iene]) ) spin[r] *= -1;
  

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
    auto spin_to_flip = (int)(uniform(gen)*(double)nspins);
    auto ix = coord[spin_to_flip][0] ;
    auto iy = coord[spin_to_flip][1] ;
    flip_spin(ix,iy);
  }

  return ;
}


//--------------------------------------------------------------------
void run_monte_carlo(int const neqsweeps, int const nsamsweeps, int const nsampstp, bool const sequential)
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
void statistics(std::vector<double> data, double &mean, double &err) {
  int ndat = data.size() ;
  double sum = 0.e0, sum2 = 0; ;
  mean = 0.e0 ;
  err=0.e0 ;

  for(int i = 0; i < ndat; i++) {
    auto val = data[i] ;
    sum += val ;
  }

  mean = sum/((double)ndat) ;

  //estimate error from jack knife data
  sum2 = 0.e0 ;
  for( int i = 0; i < ndat; i++) {
    auto delval =  (sum - data[i])/((double)(ndat-1)) - mean ;
    sum2 += delval*delval ;
  }

  
  err = sqrt(sum2*(((double)(ndat-1))/((double) ndat))) ;

  return ;
}

//--------------------------------------------------------------------
double obtain_correlations(int const tmax, std::ofstream &fp)
{
  auto ndat = observables.size() - tmax;

  if (ndat <= 0 ) {
    std::cout << observables.size() << " < tmax " << tmax << " cannot compute correlations" ;
    return 0;
  }
  
  std::vector<std::vector<double>> cor;

  std::array<double,NOBS> mean, mean2;

  for(int iobs = 0; iobs < NOBS; iobs++) {
    mean[iobs] = 0.e0;
    mean2[iobs] = 0.e0;
    for(int i = 0; i< ndat;i++) {
      auto val = observables[i][iobs] ;
      mean[iobs] += val;
      mean2[iobs] += val*val;
    }
    mean[iobs] /= (double)ndat ;
    mean2[iobs] /= (double)ndat ;
    //std::cout << iobs << " " << mean[iobs] << " " << mean2[iobs] << std::endl ;
  }


  for(int id=0; id < ndat; id++) { 
    if (id == 0) {
      for(int it=0; it < tmax; it++) {
	std::vector<double> ct;
	for(int iobs = 0; iobs < NOBS; iobs++) {
	  ct.push_back(observables[id][iobs]*observables[id+it][iobs]) ;
	}
	cor.push_back(ct);
      }
    }
    else {
      for(int it=0; it < tmax; it++) {
	for(int iobs = 0; iobs < NOBS; iobs++) {
	  cor[it][iobs] += observables[id][iobs]*observables[id+it][iobs] ;
	}
      }
    }
  }

  fp << std::scientific ; 
  for(int it=0; it < tmax; it++) {
    fp << it << " " << T << " " ;
    for(int iobs = 0; iobs < NOBS; iobs++) {
      cor[it][iobs] /= (double)ndat ;
      cor[it][iobs] = (cor[it][iobs] - mean[iobs]*mean[iobs])/
	(mean2[iobs] - mean[iobs]*mean[iobs]) ;
      fp << cor[it][iobs] << " ";
    }
    fp << std::endl ;
  }

  // Calculate the Integrated Auto-Correlation time

  std::array<double,NOBS> tau_array ;  // store autocorrelation times for each observables at a given (T, L)
  
  for(int iobs = 0; iobs < NOBS; iobs++) {
    double tau_sum = 0.50e0 ;
    for(int it=0; it < tmax; it++) {  
      tau_sum = tau_sum + cor[it][iobs];
    }
    tau_array[iobs] = tau_sum;
  }

  double tau_int = *std::max_element(tau_array.begin(), tau_array.end()); // maxima of correlation time out of all observables

  return tau_int;  // return the integrated correlation time
}


//--------------------------------------------------------------------
void analyze_samples(std::string const casename, int const nblen,
		     int const nsamsweeps, std::ofstream &fp)
{
  std::vector<std::array<double,NOBS>> meansamp, varsamp  ;
  std::array<double,NOBS> sum, sum2, avesamp, errsamp, aveflc, errflc ;
  
  auto ndat = observables.size() ;

  // std::cout << "Number of data points " << ndat << std::endl ;

  for(int iobs = 0; iobs < NOBS ; iobs++) {
    sum[iobs] = 0.e0 ;
    sum2[iobs] = 0.e0 ;
  }
  int ibin = 0;
  for(int idat = 0; idat <= ndat; idat++) {    
    if( (idat > 0 & idat%nblen == 0) || idat == ndat ) {
      for(int iobs = 0; iobs < NOBS; iobs++) {
	sum[iobs] /= (double)ibin ;
	sum2[iobs]/= (double)ibin ;
	sum2[iobs]-= (sum[iobs]*sum[iobs]) ;
	double factor;
	factor = sqrt((double(ibin))/((double)(ibin-1))) ;
	sum2[iobs] *= factor; //unbaiased estimator of variance
      }
      meansamp.push_back(sum);
      varsamp.push_back(sum2);

      //std::cout << ibin << std::endl ;
      //for(int iobs = 0; iobs < NOBS; iobs++) {
      //std::cout << iobs << " " << sum[iobs] << " " << sum2[iobs] << std::endl ;
      //}
      //std::cout << sum << std::endl ;
      //std::cout << sum2 << std::endl ;
      
      if(idat == ndat) break; 

      for(int iobs = 0; iobs < NOBS; iobs++) {
	sum[iobs] = 0.e0 ;
	sum2[iobs] = 0.e0 ;
      }
      ibin=0 ;      
    }
    for(int iobs = 0; iobs < NOBS; iobs++) {
      auto val = observables[idat][iobs] ;
      sum[iobs] += val ;
      sum2[iobs] += val*val ;
    }
    ibin++;
  }


  //std::cout << "--------" << std::endl ;
  for(int iobs = 0; iobs < NOBS; iobs++) {
    double meanval, err ;
    std::vector<double> data;
    int ndat = meansamp.size() ;
    for(int i = 0; i < ndat; i++)
      {
	data.push_back(meansamp[i][iobs]) ;
      }
    statistics(data,meanval,err) ;
    avesamp[iobs] = meanval ;
    errsamp[iobs] = err ;
    //std::cout << iobs << " " << meanval << " " << sdval << std::endl ;

    data.clear() ;
    ndat = meansamp.size() ;
    for(int i = 0; i < ndat; i++)
      {
	data.push_back(varsamp[i][iobs]) ;
      }
    statistics(data,meanval,err) ;
    aveflc[iobs] = meanval ;
    errflc[iobs] = err ;
    //std::cout << iobs << " " << meanval << " " << sdval << std::endl ;
    
  }
  

  //std::cout << "Temperature        " << T << std::endl ;
  //std::cout << "Mean Energy        " << mean[ENE] << std::endl ;
  //std::cout << "Mean Magnetization " << mean[MAG] << std::endl ;

  std::cout << T << " " << L << " ";
  std::cout << nsamsweeps << " " ;
  std::cout << avesamp[ENE] << " " << errsamp[ENE] << " " ;
  std::cout << avesamp[MAG] << " " << errsamp[MAG] << " " ;
  std::cout << aveflc[ENE]*beta*beta*double(nspins) << " " << errflc[ENE]*beta*beta*double(nspins) << " " ;
  std::cout << aveflc[MAG]*beta*(double)nspins << " " << errflc[MAG]*beta*double(nspins) << " " ;
  std::cout << avesamp[AMAG] << " " << errsamp[AMAG] << " " ;
  std::cout << aveflc[AMAG]*beta*(double)nspins << " " << errflc[AMAG]*beta*(double)nspins << " " ;
  std::cout << avesamp[MAG2] << " " << errsamp[MAG2] << " " ;
  std::cout << beta*(avesamp[MAG2] - avesamp[MAG]*avesamp[MAG])*(double)nspins  << " " << beta*(errsamp[MAG2])*(double)nspins << " " ;
  std::cout << beta*(avesamp[MAG2] - avesamp[AMAG]*avesamp[AMAG])*(double)nspins  << " " << beta*(errsamp[MAG2])*(double)nspins << " " ;
  std::cout << beta*beta*(avesamp[ENE2] - avesamp[ENE]*avesamp[ENE])*(double)nspins  << " " << beta*beta*(errsamp[ENE2])*(double)nspins << " " ;  
  std::cout << 1.e0-(avesamp[MAG4]/(3.e0*avesamp[MAG2]*avesamp[MAG2])) << " " ;
  std::cout << (avesamp[MAG2]/(avesamp[AMAG]*avesamp[AMAG])) << " " ;
  std::cout << std::endl ;

  fp << std::scientific ;
  fp << T << " " << L << " ";
  fp << nsamsweeps << " " ;
  fp << avesamp[ENE] << " " << errsamp[ENE] << " " ;
  fp << avesamp[MAG] << " " << errsamp[MAG] << " " ;
  fp << aveflc[ENE]*beta*beta*double(nspins) << " " << errflc[ENE]*beta*beta*double(nspins) << " " ;
  fp << aveflc[MAG]*beta*(double)nspins << " " << errflc[MAG]*beta*double(nspins) << " " ;
  fp << avesamp[AMAG] << " " << errsamp[AMAG] << " " ;
  fp << aveflc[AMAG]*beta*(double)nspins << " " << errflc[AMAG]*beta*(double)nspins << " " ;
  fp << avesamp[MAG2] << " " << errsamp[MAG2] << " " ;
  fp << beta*(avesamp[MAG2] - avesamp[MAG]*avesamp[MAG])*(double)nspins  << " " << beta*(errsamp[MAG2])*(double)nspins << " " ;
  fp << beta*(avesamp[MAG2] - avesamp[AMAG]*avesamp[AMAG])*(double)nspins  << " " << beta*(errsamp[MAG2])*(double)nspins << " " ;
  fp << beta*beta*(avesamp[ENE2] - avesamp[ENE]*avesamp[ENE])*(double)nspins  << " " << beta*(errsamp[ENE2])*(double)nspins << " " ;
  fp << 1.e0-(avesamp[MAG4]/(3.e0*avesamp[MAG2]*avesamp[MAG2])) << " " ;
  fp << (avesamp[MAG2]/(avesamp[AMAG]*avesamp[AMAG])) << " " ;
  fp << std::endl ;


  return ;
}


//--------------------------------------------------------------------
int main(int argc, char** argv)
{
  //PramsReader reader;
  Params parameters;
  
  std::string casename;

  casename = "tis" ;
  
  //lattice
  L=128;

  //energetics
  J=1.e0;
  T=2.2e0;
  beta=1.e0/T;

  //montecarlo
  int neqsweeps, nsamsweeps, nsampstp, nblen;
  neqsweeps = 512;
  nsamsweeps = 32768;
  nsampstp = 32 ;
  nblen = 16;

  bool sequential = true;
  
  double Tmin=2.1e0, Tmax=2.40e0, dT=0.02e0;

  bool initialize_all_up = false, print_cor = false;
  int tmax=500;
  
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
    else if( cardname == "tmax" ) {
      std::cout << "found tmax" << std::endl ;
      tmax = parameters["tmax"].asInt() ;
      std::cout << "tmax = " << tmax << std::endl ;
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

  //if(print_cor) {
    auto corfile = "output/"+casename+"_cor.plt" ;
    std::ofstream fcor(corfile) ;
    //}
  
  T = Tmin;
  
  while (T <= Tmax+dT/2.e0) { 
    //Define boltzmann weights
    define_Boltzmann();
    //initialize observable();
    initialize_observables();
    //initialize to up state();
    if(initialize_all_up) initialize_state_to_all_up() ;    

    //run
    run_monte_carlo(neqsweeps, nsamsweeps, 1, sequential) ;

    //find and print corrleations
    double tau;
    tau = obtain_correlations(tmax,fcor) ;
    std::cout << "τ = " << tau << std::endl ;

    //Define boltzmann weights
    define_Boltzmann();
    //initialize observable();
    initialize_observables();
    //initialize to up state();
    if(initialize_all_up) initialize_state_to_all_up() ;    

    //run
    run_monte_carlo(neqsweeps, nsamsweeps, 2*tau, sequential) ;

    //analyze and print output
    analyze_samples(casename, nblen, nsamsweeps,foutp) ;


    T+=dT;
  }
  foutp.close() ;
  fcor.close() ;   
    return(1) ;
}
  
