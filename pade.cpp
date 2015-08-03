//test

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <cmath>
#include <vector>
#include <complex>
#include <stdlib.h>
#include "MersenneTwister.h"
#include "fftw3.h"
using namespace std;

vector<string> split(string const& s) {
  istringstream iss(s);
  return vector<string> (istream_iterator<string>(iss),
			 istream_iterator<string>());
}

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

int main(int argc, char** argv) {
  // input from file
    ifstream fin("test.dat");
    
  char buf[1024];
  vector<vector<double> > data;
  while(fin.getline(buf,sizeof(buf))!=NULL) {
    string line = string(buf);
    vector<string> line_split = split(line);
    istringstream is1,is2;
    is1.str(line_split[0]);
    is2.str(line_split[1]);
    vector<double> line_vec(2);
    is1 >> line_vec[0];
    is2 >> line_vec[1];
    data.push_back(line_vec);
  }

  for (int i=0; i<data.size(); i++) {
    cout << data[i][0] << ' ' << data[i][1] << endl;
  }

  exit(0);

  /*
  // physical parameters
    parameters p;
    p.nsite = 256;
    p.nthermal = 5000;
    p.nbin = 100;
    p.nmeasure = 300;
    p.seed = 1231;
    p.Delta = 0.2;
    p.delta_s = 0.5;
     
    // fftw parameters
    fftw_complex *xk;
    double *xj;
    fftw_plan plan_xj_to_xk;
    
    xj=(double*)malloc(sizeof(double)*(p.nsite+2));
    xk=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(p.nsite/2+1));
    plan_xj_to_xk=fftw_plan_dft_r2c_1d(p.nsite,xj,xk,FFTW_MEASURE);
    
    fftw_execute(plan_xj_to_xk);
    // correlation func(tau)
    for(int j=0;j<p.nsite/2;j++){
       cout<<j*2*M_PI/256.0<<' '<<xk[j][0]<<endl;
    }

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
