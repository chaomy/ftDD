#ifndef _ftDD
#define _ftDD

#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <queue>
#include "armadillo"
#include "ftDDdef.hpp"

using arma::mat;
using std::cerr;
using std::copy;
using std::cout;
using std::endl;
using std::ifstream;
using std::move;
using std::ofstream;
using std::setprecision;

class ftDD {
 private:
  FitVariable dpar;  // parameters to fit
  MatConstant mcns;  // materials constant

  vector<double> ini;  // initial values
  vector<double> lob;  // lower bound of variables
  vector<double> hib;  // upper bound of variables
  vector<double> deb;  // delta = upper - lower

  vector<string> cases;  // job cases to fit

  unordered_map<string, vector<double>> weigh;  // weighs of data
  unordered_map<string, vector<double>> toltm;  // total time
  unordered_map<string, vector<double>> strss;  // total stress
  unordered_map<string, vector<double>> tstrn;  // tolal strain
  unordered_map<string, vector<double>> dnstl;  // total density
  unordered_map<string, vector<double>> flxtl;  // total flux
  unordered_map<string, vector<double>> vlc;    // total velocity
  unordered_map<string, vector<double>> ftvlc;  // fitted velocity
  unordered_map<string, DDfcc> fccslips;        // fcc density on each slips
  unordered_map<string, double> dnsprc;  // density of precipitate precinfo

  unordered_map<string, double> dparams;  // param name - params (double)
  unordered_map<string, string> sparams;  // param name - params (string)

  // fitting
  unordered_map<string,
                double (ftDD::*)(const vector<double>&, const vector<string>&)>
      calobj;
  unordered_map<string, void (ftDD::*)(const vector<string>&)> run;

 public:
  ftDD(int argc, char* argv[]);
  ~ftDD(){};

  void parseArgs(int argc, char* argv[]);
  void initParm();
  void readParam();
  void initData();
  void readData(const string&);

  // prep default setup
  void defaultSetup();

  // prep precipiates shape, distribution data
  void prepPrecDat();

  // prep DD data
  void prepData(const string&);

  // write time, stress, density, velocity
  void writeData(const string&);

  // write time, stress, density, velocity and fitted velocity
  void writeFitData(const vector<string>&);

  // write fitted variables 
  void writeResults();

  // make the fitting target more smooth (no need to use)
  void smthave(const string&);

  // make the fitting target more smooth (no need to use)
  void smthave(vector<double>&, const int);

  // load ParaDiS flux results (FCC)
  void loadDDflx(const string&);

  // load ParaDiS density (FCC)
  void loadDDdns(const string&);

  // load ParaDiS property
  void loadDDprp(const string&);

  // Using Linear mobility Law for different strian rate case
  void optErate(const vector<string>&);

  // fit all parameters, including mobility, alpha and beta
  void optPrec(const vector<string>&);

  // fix mobility, fit Alpha
  void optFixMob(const vector<string>&);

  // fix both mobility and Alpha
  void optFixMobAlpha(const vector<string>&);

  // fix all
  void optFixAll(const vector<string>&);

  // set low and high bounds for fitting variables
  void initFitting(vector<double>);

  // run CMAES algorithm do optimization
  vector<double> runCMAES(const vector<string>&);

  // calculate errors using Linear Mobility for pure metal
  double errLinearMob(const vector<double>&, const vector<string>&);

  // calculate errors using Linear Mobility for precipiates
  double errLinearPrc(const vector<double>&, const vector<string>&);

  // calculate errors using Linear Mobility use Fixed M
  double errLinearMobFixM(const vector<double>&, const vector<string>&);

  // calculate errors using Linear Mobility use Fixed M and Fixed Alpha
  double errLinearMobFixMAlp(const vector<double>&, const vector<string>&);

  // estimate fitted velocity linear mobility law
  double estLinearMob(const vector<string>&);

  void estLinearPrc(const vector<string>&);

  // utilities
  double cmaes(arma::mat& iterate, const vector<string>& kk);
  arma::mat encodev(const vector<double>& vv);
  vector<double> decodestdv(const arma::mat& vv);
};

// [a, b] -> [0, 10]
inline arma::mat ftDD::encodev(const vector<double>& vv) {
  arma::mat rs(vv.size(), 1);
  for (int i = 0; i < vv.size(); i++)
    rs[i] = 10 * acos(1. - 2. / deb[i] * (vv[i] - lob[i])) * INVPI;
  return rs;
}

// [0, 10] -> [a, b]  y = a + (b-a) × (1 – cos(π × x / 10)) / 2
inline vector<double> ftDD::decodestdv(const arma::mat& vv) {
  vector<double> rs(vv.n_elem);
  for (int i = 0; i < vv.n_elem; i++)
    rs[i] = lob[i] + deb[i] * 0.5 * (1. - cos(PI * 0.1 * vv[i]));
  return rs;
}

#endif