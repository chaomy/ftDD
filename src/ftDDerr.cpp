/*
 * @Author: chaomy
 * @Date:   2018-03-16 14:42:04
 * @Last Modified by:   chaomy
 * @Last Modified time: 2019-01-21 20:55:01
 */

#include "ftDD.hpp"

#define M 8.0

// // initialize function object
// struct errVelocityErrorFunct {
//   vector<double>& fitted_vel;  // fitted velocity
//   vector<double>& stress;      // stress vector
//   vector<double>& density;     // dislocation density
//   vector<double>& velocity;    // dislocation velocity
//   vector<double>& weight;      // weighing factors

//   errVelocityErrorFunct(vector<double>& _fitted_vel, vector<double>& _stress,
//                         vector<double>& _density, vector<double>& _velocity,
//                         vector<double>& _weight)

//       : fitted_vel(_fitted_vel),
//         stress(_stress),
//         density(_density),
//         velocity(_velocity),
//         weight(_weight) {}
// };

/*********************************************************
 * estimate fitted velocity linear mobility law
 *********************************************************/
double ftDD::estLinearMob(const vector<string>& ss) {
  double mob_burg = dpar.mob * mcns.burg;
  double alp = dpar.alp * dpar.alp, mub = mcns.mub, err = 0.0;
  for (const string& key : ss) {
    vector<double>&fv = ftvlc[key] = vector<double>(vlc[key].size());
    vector<double>&s = strss[key], &d = dnstl[key];
    vector<double>&v = vlc[key], &w = weigh[key];

    for (int i = 0; i < fv.size(); i++) {
      fv[i] = mob_burg * (s[i] - mub * sqrt(alp * d[i]));
      double r = std::abs(fv[i] - v[i]);
      err += w[i] * ((r < M) ? r * r : M * (2 * r - M));
    }
  }
  return err;
}

/*********************************************************
 * measure err using using linear mobility law
 * variables: Mob and alpha
 *********************************************************/
double ftDD::errLinearMob(const vector<double>& x, const vector<string>& ss) {
  double mob_burg = x[0] * mcns.burg;
  double alp = x[1] * x[1], mub = mcns.mub, err = 0.0;
  for (const string& key : ss) {
    vector<double>&s = strss[key], &w = weigh[key];
    vector<double>&d = dnstl[key], &v = vlc[key];
    for (int i = 0; i < s.size(); i++) {
      double r = std::abs(mob_burg * (s[i] - mub * sqrt(alp * d[i])) - v[i]);
      err += w[i] * ((r < M) ? r * r : M * (2 * r - M));
    }
  }
  return err;
}

/*********************************************************
 * estimate linear mobility law with
 *********************************************************/
void ftDD::estLinearPrc(const vector<string>& ss) {
  double mob_burg = dpar.mob * mcns.burg, mub = mcns.mub;
  double alp_alp = dpar.alp * dpar.alp;
  double bta_bta = dpar.bta * dpar.bta;

  for (const string& key : ss) {
    vector<double>&fv = ftvlc[key] = vector<double>(vlc[key].size());
    vector<double>&s = strss[key], &d = dnstl[key];
    for (int i = 0, N{static_cast<int>(s.size())}; i < N; i++)
      fv[i] = mob_burg * (s[i] - mub * sqrt(alp_alp * d[i] + bta_bta));
  }
}

/*********************************************************
 * error function, linear mobility law for precipitates
 * variables: Mob, alpha, beta
 *********************************************************/
double ftDD::errLinearPrc(const vector<double>& x, const vector<string>& ss) {
  double mob_burg = x[0] * mcns.burg;
  double alp = x[1] * x[1], mub = mcns.mub, err = 0.0;
  double dp = x[2] * x[2];
  for (const string& key : ss) {
    vector<double>&s = strss[key], &w = weigh[key];
    vector<double>&d = dnstl[key], &v = vlc[key];
    for (int i = 0; i < s.size(); i++) {
      double r =
          std::abs(mob_burg * (s[i] - mub * sqrt(alp * d[i] + dp)) - v[i]);
      err += w[i] * ((r < M) ? r * r : M * (2 * r - M));
    }
  }
  return err;
}