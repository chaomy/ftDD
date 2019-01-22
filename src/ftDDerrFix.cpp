/*
 * @Author: chaomy
 * @Date:   2018-03-16 14:42:04
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-07-11 21:34:21
 */

#include "ftDD.hpp"

#define M 8.0

/*********************************************************
 * error function, linear mobility law fixed M
 * variables: alpha
 *********************************************************/
double ftDD::errLinearMobFixM(const vector<double>& x,
                              const vector<string>& ss) {
  double mobb = dpar.mob * mcns.burg;
  double alp = x[0] * x[0], mub = mcns.mub, err = 0.0;
  for (const string& kk : ss) {
    vector<double>&s = strss[kk], &w = weigh[kk], &d = dnstl[kk], &v = vlc[kk];
    for (int i = 0; i < s.size(); i++) {
      double r = std::abs(mobb * (s[i] - mub * sqrt(alp * d[i])) - v[i]);
      err += w[i] * ((r < M) ? r * r : M * (2 * r - M));
    }
  }
  return err;
}

/*********************************************************
 * error function, linear mobility law for precipitates
 * variables: Mob, alpha, beta
 *********************************************************/
double ftDD::errLinearMobFixMAlp(const vector<double>& x,
                                 const vector<string>& ss) {
  double mobb = dpar.mob * mcns.burg;
  double alp = dpar.alp * dpar.alp;
  double mub = mcns.mub, err = 0.0;
  double dp = x[0] * x[0];
  for (const string& kk : ss) {
    vector<double>&s = strss[kk], &w = weigh[kk], &d = dnstl[kk], &v = vlc[kk];
    for (int i = 0; i < s.size(); i++) {
      double&& r = std::abs(mobb * (s[i] - mub * sqrt(alp * d[i] + dp)) - v[i]);
      err += w[i] * ((r < M) ? r * r : M * (2 * r - M));
    }
  }
  return err;
}
