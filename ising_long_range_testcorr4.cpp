// testtesttest
//test
#include "MersenneTwister.h"
#include <iostream>
#include <fstream>
#include <cmath>
// #include <vector>
using namespace std;

template<class T> inline double sqr(T x) { 
  double xtmp = (double) x;
  return xtmp*xtmp; }

template<class T> inline double quo(T x) { 
  double xtmp = (double) x;
  return xtmp*xtmp*xtmp*xtmp; }

struct parameters {
  //  int argc;
  //  char** argv;
  // parameters given by a command line
  int nsite;
  int nthermal;
  int nmeasure,nbin;
  //  int s;
  int seed;
  double Delta, delta_s;
};

class observable {
  double sum_tmp, sum_av, sum_sq;
  int ndata;
public:
  observable() {
    ndata = 0;
    sum_av = 0.0;
    sum_sq = 0.0;
  }
  ~observable() {}
  void operator<<(double x) {
    ndata++;
    sum_av += x;
    sum_sq += sqr(x);
  }
  double mean() {
    return ndata > 0 ? sum_av/(double) ndata : 0.0;
  }
  double err() {
    return ndata > 1 ? 
      sqrt((sum_sq/(double)ndata - sqr(mean()))/(double)(ndata - 1)) : 0.0;
  }
};

class Ising_MC_method {
private:
  int *spin;
  double *J,*cumulative;  
  int *cluster,*cluster_pointer;
  // for correlation function
  double *distribution;
  int *flipped_spin;
  // paremters
  parameters p;
  // random generator
  MTRand rng;
public:
  Ising_MC_method(const parameters q) :
    rng((unsigned long) q.seed), p(q) {   

    // exchange interaction
    J = new double[p.nsite];

    // probability distribution
    cumulative = new double[p.nsite];

    // spin configration    
    spin = new int[p.nsite];
    for (int i=0; i<p.nsite; i++) {
      spin[i] = 1;
    }
    // for calculation of correlation functions    
    distribution = new double[p.nsite];

    // for cluster update
    cluster.reserve(p.nsite);
    flipped_spin.reserve(p.nsite);
    // set exchange interaction (s=1 or 3)
    double J_LR, J_NN;
    // ohmic case
    J_LR = p.delta_s;
    J_NN = - log(0.5*p.Delta) - (1. + 0.577215665)*J_LR;
    J[0] = 0.0;
    for (int i=1; i<p.nsite; i++) {
      double x = M_PI*(double)i/(double)p.nsite;
      J[i] = 0.5*J_LR*V_LR(x);
    }
    J[1] += 0.5*J_NN;
    J[p.nsite-1] += 0.5*J_NN;
    // calculate cumulative distribution
    double sum_J = 0.0;
    cumulative[0] = 0.0;
    for (int i=1; i<p.nsite; i++) {
      sum_J += 2.0*J[i];
      // probability of i-th spin-flip with no spin-flip for j<i
      cumulative[i] = 1 - exp(-sum_J); 
    }
  }
  inline double V_LR(double x) {
    return sqr(M_PI/(double)p.nsite)/sqr(sin(x));
  }
  void update() {
    // update by Luijten algorithm    
    cluster.resize(0);
    flipped_spin.resize(0);

    int current_spin_initial = int(rng()*p.nsite);
    int sigma = spin[current_spin_initial];
    spin[current_spin_initial] = - spin[current_spin_initial];
    flipped_spin.push_back(current_spin_initial);
    cluster.push_back(current_spin_initial);

    while (cluster.size() > 0) {
      int current_spin = cluster.back();
      cluster.pop_back();
      
      // activate bonds
      int j=0;
      while (true) {
	// bisection search (start)
	double g = rng()*(1.0-cumulative[j]) + cumulative[j];
	if (cumulative[p.nsite-1] < g) break;
	int k_s = j+1;
	int k_e = p.nsite;
	while (k_e - k_s > 1) {
	  int k_tmp = (k_s + k_e)/2;
	  if (cumulative[k_tmp-1] < g) {
	    k_s = k_tmp;
	  } else {
	    k_e = k_tmp;
	  }
	}
	int k = k_s;
	// bisection search (end)
	j = k;
 	int i = (current_spin + k) % p.nsite;
	if (spin[i] == sigma) {
	  spin[i] = - spin[i];
	  flipped_spin.push_back(i);
	  cluster.push_back(i);
	}      
      }
    }
    // spin correlation
    for (int i=0; i<p.nsite; i++) distribution[i] = 0.0;
    for (int i=0; i<flipped_spin.size(); i++) {
      for (int j=0; j<flipped_spin.size(); j++) {
	distribution[(p.nsite + flipped_spin[i] - flipped_spin[j]) % p.nsite] += 1.;
      }
    }
    for (int i=0; i<p.nsite; i++) distribution[i] /= (double)flipped_spin.size();
  }
  double get_corr_func(const int k) {    
    return (double)distribution[k];
  }
  int get_spin(const int i) {
    return spin[i];
  }
};

int main(int argc, char** argv) {

  // physical parameters
  parameters p;
  p.nsite = 256;
  p.nthermal = 5000;
  p.nbin = 100;
  p.nmeasure = 300;
  p.seed = 1231;

  p.Delta = 0.01;
  p.delta_s = 0.5;

  // cluster update simulation
  Ising_MC_method mc(p); 
  for (int i=0; i<p.nthermal; i++) {
    mc.update();
  }    

  // observable for correlation func.
  vector<observable> corr(p.nsite,observable());
  for (int k=0; k<p.nmeasure; k++) {
    vector<observable> corr2(p.nsite,observable());
    for (int kk=0; kk<p.nbin; kk++) {
      mc.update();
      for (int j=0; j<p.nsite; j++) {
	corr2[j] << mc.get_corr_func(j);
      }
    }
    for (int j=0; j<p.nsite; j++) {
      corr[j] << corr2[j].mean();
    }
  }
  //  correlation
  for (int j=0; j<p.nsite; j++) {
    cout << j << " " 
	 << corr[j].mean() << ' ' 
	 << corr[j].err() << endl;
  }

}
