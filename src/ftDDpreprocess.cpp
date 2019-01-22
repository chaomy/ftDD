/*
 * @Author: chaomy
 * @Date:   2019-01-08 11:49:17
 * @Last Modified by:   chaomy
 * @Last Modified time: 2019-01-08 11:52:11
 */

#include "ftDD.hpp"

/*********************************************************
 * for convenience I precompute effective precipitate density
 *********************************************************/
void ftDD::prepPrecDat() {
  dnsprc["AluPrec1"] = 2.07964538e-06;
  dnsprc["AluPrec2"] = 4.12313991e-06;
  dnsprc["AluPrec3"] = 6.16509365e-06;
  dnsprc["AluPrec4"] = 8.19516017e-06;
  dnsprc["AluPrec5"] = 1.03143667e-05;

  dnsprc["AluPrecA"] = 4.39215500e-06;
  dnsprc["AluPrecB"] = 8.85722333e-06;
  dnsprc["AluPrecC"] = 1.33265420e-05;
  dnsprc["AluPrecD"] = 1.59285539e-05;
  dnsprc["AluPrecE"] = 1.77224813e-05;
}

/*********************************************************
 * make data more smooth (does nothing if width == 1)
 *********************************************************/
void ftDD::smthave(vector<double>& v, const int width) {
  double sm = 0.0;  // sliding average
  for (int i = 0; i < width; ++i) sm += v[i];
  for (int i = width; i < v.size(); i++) {
    sm += (v[i] - v[i - width]);
    v[i - width] = sm / width;
  }
}

/*********************************************************
 * add weights of each point based on the covariances
 *********************************************************/
void ftDD::smthave(const string& key) {
  int width = 1;
  smthave(vlc[key], width);
  smthave(dnstl[key], width);
  smthave(strss[key], width);

  toltm[key].erase(toltm[key].end() - width, toltm[key].end());
  vlc[key].erase(vlc[key].end() - width, vlc[key].end());
  strss[key].erase(strss[key].end() - width, strss[key].end());
  dnstl[key].erase(dnstl[key].end() - width, dnstl[key].end());
}