#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>
#include <stdlib.h>
#include "MersenneTwister.h"
#include "fftw3.h"
using namespace std;

template<class T> inline double sqr(T x) {
    double xtmp = (double) x;
    return xtmp*xtmp; }

template<class T> inline double quo(T x) {
    double xtmp = (double) x;
    return xtmp*xtmp*xtmp*xtmp; }

class cmdline {
private:
    int nint,ndouble;
    vector<int*> int_pointers;
    vector<double*> double_pointers;
    //
    string help_header;
    vector<string> type;
    vector<int> number;
    vector<string> keys,comments;
    //
    vector<string> command_line;
public:
    cmdline(const string header) : help_header(header) {
        nint = 0;
        ndouble = 0;
    }
    void add(const string key, const string comment, int &param) {
        keys.push_back(key);
        comments.push_back(comment);
        type.push_back("int");
        int_pointers.push_back(&param);
        number.push_back(nint);
        nint++;
    }
    void add(const string key, const string comment, double &param) {
        keys.push_back(key);
        comments.push_back(comment);
        type.push_back("double");
        double_pointers.push_back(&param);
        number.push_back(ndouble);
        ndouble++;
    }
    void help() {
        cout << help_header << endl;
        for (int i=0; i<keys.size(); i++) {
            cout << keys[i] << " (" << type[i] <<"): " << comments[i] << endl;
        }
    }
    void parse(const int argc, char** argv) {
        if (argc <= 1) {
            help();
            exit(1);
        }
        for (int i=1; i<=argc-1; i++) {
            command_line.push_back(argv[i]);
        }
        for (int i=0; i<keys.size(); i++) {
            for (int j=0; j<command_line.size(); j++) {
                if (keys[i] == command_line[j]) {
                    if (type[i] == "int") {
                        *int_pointers[number[i]] = atoi(command_line[j+1].c_str());
                    } else if (type[i] == "double") {
                        *double_pointers[number[i]] = atof(command_line[j+1].c_str());
                    } 
                    break;
                } 
                if (j == command_line.size() - 1) {
                    help();
                    exit(1);
                }
            }
        }
    }
};

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
    vector<int> spin,cluster;
    vector<double> J,cumulative;
    // for correlation function
    vector<int> flipped_spin;
    vector<double> distribution;
    // paremters
    parameters p;
    // random generator
    MTRand rng;
public:
    Ising_MC_method(const parameters q) :
    rng((unsigned long) q.seed), p(q) {
        
        // exchange interaction
        J.reserve(p.nsite);
        // probability distribution
        cumulative.reserve(p.nsite);
        // spin configration
        spin.reserve(p.nsite);
        for (int i=0; i<p.nsite; i++) {
            spin[i] = 1;
        }
        // for calculation of correlation functions
        distribution.reserve(p.nsite);
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
    bool command_line = true;
    
    if (command_line) {
        cmdline cmd("Options");
        cmd.add("nthermal","# of MCS for thermalization",p.nthermal);
        cmd.add("nbin","# of bins",p.nbin);
        cmd.add("nmeasure","# of measurement",p.nmeasure);
        // MCS = nbin*nmeasure
        cmd.add("delta","tunneling matrix",p.Delta);
        cmd.add("nsite","# of sites (=beta)",p.nsite);
        cmd.add("alpha","dissipation strength",p.delta_s);
        cmd.add("nseed","seed of random number generator",p.seed);
        cmd.parse(argc,argv);
    } else {
        p.nsite = 256;
        p.nthermal = 5000;
        p.nbin = 100;
        p.nmeasure = 300;
        p.seed = 1231;
        p.Delta = 0.2;
        p.delta_s = 0.5;
    }

    /*
    cout << p.nsite << ' '
         << p.Delta << ' '
         << p.seed << ' '
        << p.nsite << ' ' << p.delta_s << ' '
        << p.nthermal << ' ' << p.nbin << ' ' << p.nmeasure << endl;
    */
    //    exit(0);
    
     
    // fftw parameters
    /*
    fftw_complex *xk;
    double *xj;
    fftw_plan plan_xj_to_xk;
    
    xj=(double*)malloc(sizeof(double)*(p.nsite+2));
    xk=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(p.nsite/2+1));
    plan_xj_to_xk=fftw_plan_dft_r2c_1d(p.nsite,xj,xk,FFTW_MEASURE);
    */    

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
    
    // correlation
    for (int j=0; j<p.nsite; j++) {
      cout << (double)j << ' '
	   << corr[j].mean() << ' '
	   << corr[j].err() << endl;
    }

    /*
    
    fftw_execute(plan_xj_to_xk);
    // correlation func(tau)
    for(int j=0;j<p.nsite/2;j++){
       cout<<j*2*M_PI/256.0<<' '<<xk[j][0]<<endl;
    }
    exit(0);

    const int ndata = p.nsite/4; // ndata must be larger than 2
   
    // Pade approximation
    // Pade parameters
    std::complex<double> a[ndata],z[ndata],**g,A[ndata+1],B[ndata+1];
    
    // Matsubara frequency on Im axis
    for(int i=0;i<ndata;i++){
        z[i]=std::complex<double> (0,2.0*M_PI/p.nsite*i);
    }
    
    g=new complex<double>* [ndata];
    for(int i=0;i<ndata;i++){
        g[i]=new complex<double>[i+1];
    }
    
    // set the initial value of reccurence relation
    for(int i=0;i<ndata;i++){
        g[i][0]=xk[i][0];
    }
    
    // calculate reccurence relation
    for(int i=1;i<ndata;i++){
        for(int j=1;j<=i;j++){
            g[i][j]=(g[j-1][j-1]/g[i][j-1]-std::complex<double>(1,0))/(z[i]-z[j-1]);
        }
    }
    for(int i=0;i<ndata;i++){
        a[i]=g[i][i];
    }
    
    double coeff = M_PI/4.*p.Delta*p.Delta;
    // evaluate spectrum func.
    double xs = 0.;
    double xe = 5.;
    int nd = 200;
    double xh = (xe - xs)/(double) nd;
    
    for (int j=0; j<=nd; j++) {
      double x = xs + xh*(double) j;
      complex<double> tmp(1.,0.);
      for (int p=ndata-1; p>0; p--) {
          tmp = 1.0 + a[p]*(complex<double>(x*coeff,0.)-z[p-1])/tmp;
      }
      tmp = a[0]/tmp;
      cout << x << ' '  << tmp.imag()/x/43 <<

    fftw_destroy_plan(plan_xj_to_xk);
    fftw_free(xk);fftw_free(xj);
    */
    
}
