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
  FitVariable dpar;  // fitted parameters
  MatConstant mcns;  // materials constant

  vector<double> ini, lob, hib, deb;
  vector<string> cases;  // cases to fit

  unordered_map<string, vector<double>> weigh;  // weighs of data
  unordered_map<string, vector<double>> toltm;  // total time
  unordered_map<string, vector<double>> strss;  // total stress
  unordered_map<string, vector<double>> tstrn;  // tolal strain
  unordered_map<string, vector<double>> dnstl;  // total density
  unordered_map<string, vector<double>> flxtl;  // total flux
  unordered_map<string, vector<double>> vlc;    // total velocity
  unordered_map<string, vector<double>> ftvlc;  // fitted velocity
  unordered_map<string, DDfcc> dd;
  unordered_map<string, double> dnsprc;  // map type to precinfo

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

  // write fitting variables
  void writeResults();

  // to make the fitting target more smooth (not used)
  void smthave(const string&);

  // to make the fitting target more smooth (not used)
  void smthave(vector<double>&, const int);

  // load ParaDiS flux results (FCC)
  void loadDDflx(const string&);

  // load ParaDiS density (FCC)
  void loadDDdns(const string&);

  // load ParaDiS property
  void loadDDprp(const string&);

  // Using Linear mobility Law for different strian rate case
  void optErate(const vector<string>&);

  //
  void optPrec(const vector<string>&);

  //
  void optFixMob(const vector<string>&);

  //
  void optFixMobAlpha(const vector<string>&); 

  //
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