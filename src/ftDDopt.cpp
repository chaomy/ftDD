/*
 * @Author: chaomy
 * @Date:   2018-03-16 14:42:04
 * @Last Modified by:   chaomy
 * @Last Modified time: 2019-01-08 12:33:13
 */

#include "ftDD.hpp"

/*********************************************************
 * set low and high bounds for fitting variables
 *********************************************************/
void ftDD::initFitting(vector<double> v) {
  ini = v, lob = v, hib = v, deb = v;
  for (int i = 0; i < ini.size(); ++i)
    deb[i] = (hib[i] *= 10000) - (lob[i] /= 10000);
}

/*********************************************************
 * run CMAES algorithm to optimize variables
 *********************************************************/
vector<double> ftDD::runCMAES(const vector<string>& keys) {
  arma::mat it(encodev(ini));
  dpar.err = cmaes(it, keys);
  return decodestdv(it);
}

/*********************************************************
 * fix mob and do the fitting on alpha
 *********************************************************/
void ftDD::optFixMob(const vector<string>& keys) {
  initFitting({dpar.alp});
  vector<double> res(move(runCMAES(keys)));
  dpar.alp = res[0];
  dpar.err = estLinearMob(keys);
  writeFitData(keys);
}

/*********************************************************
 * fix Mob and alpha, and do the fitting on beta
 *********************************************************/
void ftDD::optFixMobAlpha(const vector<string>& keys) {
  initFitting({dpar.bta});
  vector<double> res(move(runCMAES(keys)));
  dpar.bta = res[0];
  estLinearPrc(keys);
  writeFitData(keys);
}

/*********************************************************
 * fix all, used to check errors of results
 *********************************************************/
void ftDD::optFixAll(const vector<string>& keys) {
  dpar.err = estLinearMob(keys);
  writeFitData(keys);
}

/*********************************************************
 * Using Linear mobility Law for different strian rate case
 * fitting variables: mob and alpha
 *********************************************************/
void ftDD::optErate(const vector<string>& keys) {
  // initialize fitting parameters
  initFitting({dpar.mob, dpar.alp});
  vector<double> res(std::move(runCMAES(keys)));
  dpar.mob = res[0], dpar.alp = res[1];
  dpar.err = estLinearMob(keys);
  writeFitData(keys);
}

/*********************************************************
 * fitting on precipitates
 *********************************************************/
void ftDD::optPrec(const vector<string>& keys) {  
  initFitting({dpar.mob, dpar.alp, dpar.bta});
  vector<double> res(move(runCMAES(keys)));
  dpar.mob = res[0], dpar.alp = res[1], dpar.bta = res[2];
  estLinearPrc(keys);
  writeFitData(keys);
}